#include "stdafx.h"
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
// Support vector machine optimisation functions
//
// Written by: Alistair Shilton
//             Melbourne University
//

#include "svoptim.h"
#include "svdata.h"
#include "svflags.h"
#include "common/vector.h"
#include "common/matrix.h"
#include "common/factor.h"

unsigned long solve_active(SVdata &problem, L_DOUBLE &sol_tol, unsigned long epochs, volatile int *async_exit_flag, volatile unsigned long *loop_count);
unsigned long solve_active_fixed_bias(SVdata &problem, L_DOUBLE &sol_tol, unsigned long epochs, volatile int *async_exit_flag, volatile unsigned long *loop_count);

void calc_step_active(SVdata &problem, long &n_zero, int f_zero, fVECTOR &d_alpha_pivot, L_DOUBLE &d_b, fVECTOR &d_e_pivot, L_DOUBLE &d_f, long &start, long &finish);
void calc_step_scale_active(SVdata &problem, long &p, int &p_indic, fVECTOR &d_alpha_pivot, L_DOUBLE &d_b, fVECTOR &d_e_pivot, L_DOUBLE &d_f, long start, long finish);
void take_step_active(SVdata &problem, long &n_zero, int &f_zero, int p_indic, fVECTOR &d_alpha_pivot, L_DOUBLE &d_b, fVECTOR &d_e_pivot, L_DOUBLE &d_f, long start, long finish);
int test_optimality_active(SVdata &problem, long &q, int &q_indic, int f_zero, int p_indic, L_DOUBLE &sol_tol);
void fix_active_set_active(SVdata &problem, long p, int p_indic, long q, int q_indic, long &n_zero);
void calc_step_active_fixed_bias(SVdata &problem, long &n_zero, fVECTOR &d_alpha_pivot, fVECTOR &d_e_pivot, long &start, long &finish);
void calc_step_scale_active_fixed_bias(SVdata &problem, long &p, int &p_indic, fVECTOR &d_alpha_pivot, fVECTOR &d_e_pivot, long start, long finish);
void take_step_active_fixed_bias(SVdata &problem, long &n_zero, int p_indic, fVECTOR &d_alpha_pivot, fVECTOR &d_e_pivot, long start, long finish);
int test_optimality_active_fixed_bias(SVdata &problem, long &q, int &q_indic, int p_indic, L_DOUBLE &sol_tol);
void fix_active_set_active_fixed_bias(SVdata &problem, long p, int p_indic, long q, int q_indic, long &n_zero);

unsigned long solve_SMO(SVdata &problem, L_DOUBLE &sol_tol, unsigned long epochs, volatile int *async_exit_flag, volatile unsigned long *loop_count, JTYPE *J);
unsigned long solve_SMO_fixed_bias(SVdata &problem, L_DOUBLE &sol_tol, unsigned long epochs, volatile int *async_exit_flag, volatile unsigned long *loop_count);

int examineExample_SMO(long i2, int tau2, int f_zero, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J);
long secondChoice_SMO(L_DOUBLE &E1, L_DOUBLE &E2, SVdata &problem);
int takeStep_SMO(long i1, long i2, int tau1, int tau2, L_DOUBLE &E1, L_DOUBLE &E2, int f_zero, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J);
int trial_step_SMO(L_DOUBLE &d_J_epart, long i1, long i2, int tau1, int tau2, L_DOUBLE &e1, L_DOUBLE &e2, L_DOUBLE &d_alpha1, L_DOUBLE &d_alpha2, L_DOUBLE &d_b, L_DOUBLE &d_f, int &f_zero, int &f_zero_next, SVdata &problem, L_DOUBLE &sol_tol);
int actually_take_step_SMO(L_DOUBLE &d_J_epart, long i1, long i2, int tau1, int tau2, L_DOUBLE &d_alpha1, L_DOUBLE &d_alpha2, L_DOUBLE &d_b, L_DOUBLE &d_f, int &f_zero, int &f_zero_next, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J);
int takeStep_SMO_fixed_bias(long i, int tau, SVdata &problem, L_DOUBLE &sol_tol);

unsigned long solve_d2c(SVdata &problem, L_DOUBLE &sol_tol, unsigned long epochs, volatile int *async_exit_flag, volatile unsigned long *loop_count, JTYPE *J);
unsigned long solve_d2c_fixed_bias(SVdata &problem, L_DOUBLE &sol_tol, unsigned long epochs, volatile int *async_exit_flag, volatile unsigned long *loop_count);

int takeStep_d2c(long i1, long i2, int tau1, int tau2, L_DOUBLE &e1, L_DOUBLE &e2, int f_zero, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J);
int takeStep_d2c_pattern(long i1, long i2, int tau1, int tau2, L_DOUBLE &e1, L_DOUBLE &e2, int f_zero, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J);
int takeStep_d2c_pattern_fastkern(long i1, long i2, int tau1, int tau2, L_DOUBLE &e1, L_DOUBLE &e2, int f_zero, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J);
int trial_step_d2c(L_DOUBLE &d_J_epart, long i1, long i2, int tau1, int tau2, L_DOUBLE &e1, L_DOUBLE &e2, L_DOUBLE &d_alpha1, L_DOUBLE &d_alpha2, L_DOUBLE &d_b, L_DOUBLE &d_f, int &f_zero, int &f_zero_next, SVdata &problem, L_DOUBLE &sol_tol);
int trial_step_d2c_fastkern(L_DOUBLE &d_J_epart, long i1, long i2, int tau1, int tau2, L_DOUBLE &e1, L_DOUBLE &e2, L_DOUBLE &d_alpha1, L_DOUBLE &d_alpha2, L_DOUBLE &d_b, L_DOUBLE &d_f, int &f_zero, int &f_zero_next, SVdata &problem, L_DOUBLE &sol_tol);
int actually_take_step_d2c(L_DOUBLE &d_J_epart, long i1, long i2, int tau1, int tau2, L_DOUBLE &d_alpha1, L_DOUBLE &d_alpha2, L_DOUBLE &d_b, L_DOUBLE &d_f, int &f_zero, int &f_zero_next, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J);
int actually_take_step_d2c_pattern(L_DOUBLE &d_J_epart, long i1, long i2, int tau1, int tau2, L_DOUBLE &d_alpha1, L_DOUBLE &d_alpha2, L_DOUBLE &d_b, L_DOUBLE &d_f, int &f_zero, int &f_zero_next, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J);



#define E_CHECKER                                                           \
{                                                                           \
        L_DOUBLE __new_e;                                                   \
        long __I__pivot;                                                    \
        long __J__pivot;                                                    \
        std::cerr << "phantom - echeck\n";                                  \
        for ( __I__pivot = 1 ; __I__pivot <= problem.N ; __I__pivot++ )     \
        {                                                                   \
            __new_e = problem.b - (problem.fast_z)[(problem.fast_m)[__I__pivot-1]-1]; \
                                                                            \
            if ( problem.N_Z < problem.N )                                  \
            {                                                               \
                for ( __J__pivot = (problem.N_Z)+1 ; __J__pivot <= problem.N ; __J__pivot++ ) \
                {                                                           \
                    __new_e += (problem.opt_get_G_tau_pivot(__I__pivot-1,__J__pivot-1)) * (problem.fast_alpha)[(problem.fast_m)[__J__pivot-1]-1]; \
                }                                                           \
            }                                                               \
                                                                            \
            if ( (problem.fast_tau)[(problem.fast_m)[__I__pivot-1]-1] > 0 )           \
            {                                                               \
                __new_e += (problem.fast_rho)[(problem.fast_m)[__I__pivot-1]-1];      \
            }                                                               \
                                                                            \
            else if ( (problem.fast_tau)[(problem.fast_m)[__I__pivot-1]-1] < 0 )      \
            {                                                               \
                __new_e -= (problem.fast_rho_star)[(problem.fast_m)[__I__pivot-1]-1]; \
            }                                                               \
                                                                            \
            std::cerr << __I__pivot << ": " << problem.opt_get_e_pivot(__I__pivot-1) << " - " << __new_e << " = " << problem.opt_get_e_pivot(__I__pivot-1) - __new_e << "\n"; \
                                                                            \
            if ( fabs((problem.opt_get_e_pivot(__I__pivot-1))-__new_e) > 1e-6 )   \
            {                                                               \
                std::cerr << "fuck fuck fuck fuck\n";                       \
                exit(1);                                                    \
            }                                                               \
        }                                                                   \
}

#define E_CHECKER_QUIET                                                     \
{                                                                           \
        L_DOUBLE __wen_f_;                                                  \
        long __I__pviot_;                                                   \
        long __J__pviot_;                                                   \
        std::cerr << "phantom - echeck\n";                                  \
        for ( __I__pviot_ = 1 ; __I__pviot_ <= problem.N ; __I__pviot_++ )  \
        {                                                                   \
            __wen_f_ = problem.b - (problem.fast_z)[(problem.fast_m)[__I__pviot_-1]-1]; \
                                                                            \
            if ( problem.N_Z < problem.N )                                  \
            {                                                               \
                for ( __J__pviot_ = (problem.N_Z)+1 ; __J__pviot_ <= problem.N ; __J__pviot_++ ) \
                {                                                           \
                    __wen_f_ += (problem.opt_get_G_tau_pivot(__I__pviot_-1,__J__pviot_-1)) * (problem.fast_alpha)[(problem.fast_m)[__J__pviot_-1]-1]; \
                }                                                           \
            }                                                               \
                                                                            \
            if ( (problem.fast_tau)[(problem.fast_m)[__I__pviot_-1]-1] > 0 )          \
            {                                                               \
                __wen_f_ += (problem.fast_rho)[(problem.fast_m)[__I__pviot_-1]-1];    \
            }                                                               \
                                                                            \
            else if ( (problem.fast_tau)[(problem.fast_m)[__I__pviot_-1]-1] < 0 )     \
            {                                                               \
                __wen_f_ -= (problem.fast_rho_star)[(problem.fast_m)[__I__pviot_-1]-1]; \
            }                                                               \
                                                                            \
            if ( fabs((problem.opt_get_e_pivot(__I__pviot_-1))-__wen_f_) > 1e-6 ) \
            {                                                               \
                std::cerr << "fuck fuck fuck fuck\n";                       \
                E_CHECKER;                                                  \
                exit(1);                                                    \
            }                                                               \
        }                                                                   \
}


unsigned long solve(SVdata &problem, L_DOUBLE &sol_tol, unsigned long epochs, volatile int *async_exit_flag, volatile unsigned long *loop_count, volatile unsigned long *loop_count_last, JTYPE *J)
{
    long i;
    long n_neg;
    long n_pos;
    unsigned long result = 0;
    L_DOUBLE d_b;
    L_DOUBLE temp;

    (*loop_count) = 1;

    switch ( problem.N )
    {
        case 0:
        {
            // The problem is ill-defined, so just leave the bias b
            // as whatever it is now and get out.

            #ifdef DEBUG
            std::cerr << "problem size 0\n";
            #endif

            break;
        }

        default:
        {
            problem.enter_opt();

            // We need to ensure that the problem is well defined in a
            // geometric sense.  That is, if all our training points lie on
            // the same side of the decision surface then the "solution" is
            // just to choose an arbitrary w (equivalently alpha) (we use
            // whatever it is now) and then let b grow to +/- infinity
            // (remember that b is NOT bounded in the SVM optimisation
            // problem).  To overcome this, we just choose an appropriate
            // value for b to satisfy all the constraints, given the current
            // value of w (alpha).
            //
            // NB: - if the bias is fixed, this won't be a problem, so we
            //       can just enter the problem as usual.
            //     - if *all* alpha's are constrained to 0 then there is
            //       no problem (as there are no constraints on e in this
            //       case - e can be *anything* and it's "optimal" - so
            //       we will just stick with whatever bias b we have).

            n_neg = 0;
            n_pos = 0;

            #ifdef DEBUG
            std::cerr << "problem size " << problem.N << "\n";
            #endif

            for ( i = 1 ; ( ( i <= problem.N ) && ( ( n_neg == 0 ) || ( n_pos == 0 ) ) ) ; i++ )
            {
                switch ( (problem.fast_contype)[i-1] )
                {
                    case 0:  {                   break; }
                    case 1:  { n_neg++;          break; }
                    case 2:  {          n_pos++; break; }
                    default: { n_neg++; n_pos++; break; }
                }
            }

            // If n_neg = n_pos = 0 then all alpha is constrained to zero,
            // so there really is nothing to be done here (except possibly
            // to set b appropriately, but how are you going to do that
            // for such an ill-defined problem?).

            if ( ( n_neg != 0 ) || ( n_pos != 0 ) )
            {
                // OK, we have a problem.  The question is - is it well
                // defined?  If it is ill-defined and fixed bias, then
                // there is nothing much to be done.  However, the standard
                // solver can deal with it (in the variable bias case, this
                // will just send the bias to +-infinity, which is not a
                // good thing).

                #ifdef DEBUG
                std::cerr << "Have points\n";
                #endif

                if ( ( ( n_neg == 0 ) || ( n_pos == 0 ) ) && !(problem.fix_bias) )
                {
                    // Problem is ill-defined, variable bias.  Some alpha's
                    // may be non-zero, some may not.  The current bias may
                    // be non-zero, or it may not.  This doesn't matter.
                    // We just need to find a bias step that is sufficient
                    // to ensure that e is correct.
                    //
                    // all +: get_e_pos(i) >= 0
                    // all -: get_e_neg(i) <= 0
                    //
                    // where e_i is what is reported by the SVdata problem.
                    // Hence all the training points will be correctly
                    // classified, whatever that means.  This code deals
                    // with this by finding the minimum b change to make
                    // this true.

                    #ifdef DEBUG
                    std::cerr << "Problem ill-defined.  Setting b... ";
                    std::cerr << problem << "\n\n";
                    #endif

                    d_b = 0.0;

                    if ( n_neg == 0 )
                    {
                        for ( i = 1 ; i <= problem.N ; i++ )
                        {
                            temp = problem.opt_get_e_pos_pivot(i-1);

                            if ( -temp >= d_b )
                            {
                                d_b = -temp;
                            }
                        }
                    }

                    else
                    {
                        for ( i = 1 ; i <= problem.N ; i++ )
                        {
                            temp = problem.opt_get_e_neg_pivot(i-1);

                            if ( -temp <= d_b )
                            {
                                d_b = -temp;
                            }
                        }
                    }

                    problem.opt_step_none(d_b);

                    #ifdef DEBUG
                    std::cerr << "Finished ill-defined steppage.\n";
                    std::cerr << problem << "\n\n";
                    #endif
                }

                else
                {
                    // The problem is well defined, so we can palm it off to
                    // our optimiser of choice.

                    #ifdef DO__FLOPS
                    SET_FLOP_MODE_OPTIM;
                    #endif

                    #ifdef DEBUG
                    std::cerr << "Jump to optimiser\n";
                    #endif

                    if ( problem.fix_bias )
                    {
                        if ( OPT_ACTIVE_SET(problem.svflags) )
                        {
                            #ifdef DEBUG
                            std::cerr << "Fixed bias active\n";
                            #endif

                            result = solve_active_fixed_bias(problem,sol_tol,epochs,async_exit_flag,loop_count);
                        }

                        else if ( OPT_PLATT_SMO(problem.svflags) )
                        {
                            default_fixbias_opt:

                            #ifdef DEBUG
                            std::cerr << "Fixed bias smo\n";
                            #endif

                            result = solve_SMO_fixed_bias(problem,sol_tol,epochs,async_exit_flag,loop_count);
                        }

                        else if ( OPT_DANIEL_D2C(problem.svflags) )
                        {
                            #ifdef DEBUG
                            std::cerr << "Fixed bias d2c\n";
                            #endif

                            result = solve_d2c_fixed_bias(problem,sol_tol,epochs,async_exit_flag,loop_count);
                        }

                        else
                        {
                            goto default_fixbias_opt;
                        }
                    }

                    else
                    {
                        if ( OPT_ACTIVE_SET(problem.svflags) )
                        {
                            #ifdef DEBUG
                            std::cerr << "variable bias active\n";
                            #endif

                            result = solve_active(problem,sol_tol,epochs,async_exit_flag,loop_count);
                        }

                        else if ( OPT_PLATT_SMO(problem.svflags) )
                        {
                            default_varbias_opt:

                            #ifdef DEBUG
                            std::cerr << "variable bias smo\n";
                            #endif

                            result = solve_SMO(problem,sol_tol,epochs,async_exit_flag,loop_count,J);
                        }

                        else if ( OPT_DANIEL_D2C(problem.svflags) )
                        {
                            #ifdef DEBUG
                            std::cerr << "variable bias d2c\n";
                            #endif

                            result = solve_d2c(problem,sol_tol,epochs,async_exit_flag,loop_count,J);
                        }

                        else
                        {
                            goto default_varbias_opt;
                        }
                    }

                    #ifdef DO__FLOPS
                    SET_FLOP_MODE_MISC;
                    #endif
                }
            }

            problem.exit_opt();

            break;
        }
    }

    (*loop_count_last) += (*loop_count);
    (*loop_count)       = 0;

    return result;
}



/***********************************************************************

                        Active set optimiser
                        ====================

***********************************************************************/


unsigned long solve_active(SVdata &problem, L_DOUBLE &sol_tol, unsigned long epochs, volatile int *async_exit_flag, volatile unsigned long *loop_count)
{
    long i;

    #ifdef DEBUG
    long i_pivot;
    #endif

    long start;
    long finish;

    long n_zero;
    int f_zero;

    fVECTOR d_alpha_pivot('r',problem.N);
    fVECTOR d_e_pivot('r',problem.N);
    L_DOUBLE d_b;
    L_DOUBLE d_f;

    long p = 0;
    int p_indic = 0;

    long q = 0;
    int q_indic = 0;

    //
    // n_zero = number of zero's starting e_b
    //
    // e_b = [ e_bz  ] = [  0    ]  (vector containing n_zero elements)
    //       [ e_bnz ]   [ e_bnz ]  (vector containing N-n_zero elements)
    //
    // f_zero = 0 => f != 0, so we still need to deal with this.
    //        = 1 => f == 0, and hence will be for the rest of the alg.
    //
    // Scaling factor and marker:
    //
    // p       = variable that has hit bound to cause
    //           scaling (0 if none).
    // p_indic = 0  if no variable has hit bound.
    //           +1 hit 0 bound from above.
    //           -1 hit v bound from above.
    //           +2 hit h bound from below.
    //           -2 hit 0 bound from below.
    //
    // Least non-optimal gradient and marker of variable:
    //
    // q       = variable that is most non-optimal (0=none).
    // q_indic = 0  if all constrained variables are optimal.
    //           +1 variable (at 0) will improve increasing.
    //           +2 variable (at h) will improve decreasing.
    //           -1 variable (at 0) will improve decreasing.
    //           -2 variable (at v) will improve increasing.
    //

    //
    // Set initial state for f_zero and n_zero
    //

    if ( problem.f == 0.0 ) { f_zero = 1; }
    else                    { f_zero = 0; }

    n_zero = 0;

    if ( problem.N_F > 0 )
    {
        for ( i = (problem.N_C)+1 ; i <= problem.N ; i++ )
        {
            if ( problem.opt_get_e_pivot(i-1) != 0.0 )
            {
                break;
            }

            n_zero++;
        }
    }

    //
    // Main algorithmic loop.
    //

    for ( ; ( ( (*loop_count) <= epochs ) && ( !(*async_exit_flag) ) ) ; (*loop_count)++ )
    {
        #ifdef DEBUG
        E_CHECKER
        #endif   

        #ifdef DEBUG
        std::cerr << "===========" << (*loop_count) << "=============\n";
        {
            SVdata_short_summary temprob(problem);
            std::cerr << temprob << "\n\n";
        }
        #endif

        #ifdef DO__FLOPS
        REG_n_iterations;
        #endif

        //
        // Calculate the step (ignoring bounds on free variables).  This
        // will also detect (and fix) any changes to n_zero.
        //

        start  = (problem.N_C)+1;
        finish = problem.N;

        #ifdef DEBUG
        std::cerr << "phantom step a\n";
/*phantom*/
std::cerr << problem.G_tau << "\n";
        #endif

        calc_step_active(problem,n_zero,f_zero,d_alpha_pivot,d_b,d_e_pivot,d_f,start,finish);

        //
        // Calculate scale factor for step.  The scale factor is deta, and
        // p and p_indic mark the alpha which caused this as described
        // previously.
        //

        #ifdef DEBUG
        std::cerr << "Phantom step b\n";
        std::cerr << "d_alpha: " << d_alpha_pivot << "\n\n";
        std::cerr << "d_e: " << d_e_pivot << "\n\n";
        std::cerr << "d_b: " << d_b << "\n\n";
        std::cerr << "d_f: " << d_f << "\n\n";
        if ( problem.N_F )
        {
            long iq,ir;
            std::cerr << "G_tau: \n\n";
            for ( iq = ((problem.N_C)+1) ; iq <= problem.N ; iq++ )
            {
                for ( ir = ((problem.N_C)+1) ; ir <= problem.N ; ir++ )
                {
                    std::cerr << (problem.opt_get_G_tau_pivot(iq-1,ir-1)) << "\t";
                }
                std::cerr << "\n";
            }
            std::cerr << "\n";
        }
        #endif

        calc_step_scale_active(problem,p,p_indic,d_alpha_pivot,d_b,d_e_pivot,d_f,start,finish);

        //
        // Take scaled step, fix f_zero and n_zero appropriately.
        //

        #ifdef DEBUG
        std::cerr << "Phantom step c\n";
        std::cerr << "d_alpha: " << d_alpha_pivot << "\n\n";
        std::cerr << "d_e: " << d_e_pivot << "\n\n";
        std::cerr << "d_b: " << d_b << "\n\n";
        std::cerr << "d_f: " << d_f << "\n\n";
        #endif

        take_step_active(problem,n_zero,f_zero,p_indic,d_alpha_pivot,d_b,d_e_pivot,d_f,start,finish);

        #ifdef DEBUG
        std::cerr << "===========" << (*loop_count) << " phantom=============\n";
        std::cerr << "d_alpha: " << d_alpha_pivot << "\n\n";
        std::cerr << "d_e: " << d_e_pivot << "\n\n";
        std::cerr << "d_b: " << d_b << "\n\n";
        std::cerr << "d_f: " << d_f << "\n\n";
        std::cerr << "start = " << start << "\n";
        std::cerr << "finish = " << finish << "\n";
        {
            SVdata_short_summary temprob(problem);
            std::cerr << temprob << "\n\n";
        }
        #endif

        //
        // Check optimality of solution, exit if optimal.  Also find the
        // most inoptimal variable if solution is not optimal.
        //

        #ifdef DEBUG
        std::cerr << "Phantom step d\n";
        #endif

        if ( test_optimality_active(problem,q,q_indic,f_zero,p_indic,sol_tol) )
        {
            break;
        }

        //
        // Modify the active set.  Fix n_zero if required.
        //

        #ifdef DEBUG
        std::cerr << "Phantom step e\n";
        #endif

        fix_active_set_active(problem,p,p_indic,q,q_indic,n_zero);
    }

    return (*loop_count);
}

void calc_step_active(SVdata &problem, long &n_zero, int f_zero, fVECTOR &d_alpha_pivot, L_DOUBLE &d_b, fVECTOR &d_e_pivot, L_DOUBLE &d_f, long &start, long &finish)
{
    long nbad;

    L_DOUBLE theta;
    L_DOUBLE temp;

    //
    // Default step is zero
    //

    #ifdef DEBUG
    std::cerr << "zeroing step to begin\n";
    #endif

    d_b = 0.0;
    d_f = 0.0;

    d_e_pivot     = 0.0;
    d_alpha_pivot = 0.0;

    //
    // We only calculate a step if N_F > 0
    //

    if ( problem.N_F > 0 )
    {
        //
        // Find how much (if any) of V is singular
        //
        // Note that if fixed bias is selected then V functions will ensure
        // that d_b is set appropriately to zero.
        //

        #ifdef DEBUG
        std::cerr << "size nonzero\n";
        #endif

        if ( ( nbad = problem.opt_get_nbad_factor() ) == 0 )
        {
            //
            // V is non-singular, so step is Newton.
            //

            #ifdef DEBUG
            std::cerr << "starting nonsing (trivial?) step\n";
            #endif

            start  = (problem.N_C)+1;
            finish = problem.N;

            #ifdef DO__FLOPS
            REG_n_step_nonsing;
            #endif

            //
            // If e is all zero, the step is trivial so no need to bother
            // actually calculating it.
            //

            if ( n_zero < problem.N_F )
            {
                #ifdef DEBUG
                std::cerr << "starting nonsing nontriv step\n";
                #endif

                d_e_pivot.set_offset_start((problem.N_C)+1);
                d_e_pivot.set_offset_end(problem.N);
                d_alpha_pivot.set_offset_start((problem.N_C)+1);
                d_alpha_pivot.set_offset_end(problem.N);

                #ifdef DEBUG
                std::cerr << "minverse step is here\n";
                #endif

                problem.opt_minverse_factor(d_alpha_pivot,d_b,n_zero,0,f_zero);

                #ifdef DEBUG
                std::cerr << "free pivot alloc\n";
                #endif

                d_e_pivot = problem.opt_get_e_free_pivot();
                d_f = problem.f;

                #ifdef DEBUG
                std::cerr << "negation point\n";
                #endif

                d_alpha_pivot.negate_it();
                d_e_pivot.negate_it();

                d_b = -d_b;
                d_f = -d_f;

                d_alpha_pivot.reset_offsets();
                d_e_pivot.reset_offsets();

                #ifdef DEBUG
                std::cerr << "phantom phantom\n";
                std::cerr << start << " - " << finish << "\n";
                std::cerr << (problem.N_C) << " - " << problem.N << "\n";
                std::cerr << "d_e: " << d_e_pivot << "\n\n";
                std::cerr << "d_alpha: " << d_alpha_pivot << "\n\n";
                std::cerr << "d_b: " << d_b << "\n\n";
                std::cerr << "d_f: " << d_f << "\n\n";
                #endif
            }
        }
    
        else
        {
            //
            // V is singular, so we step in a direction of linear non-ascent
            // w.r.t. alpha (ignoring possible non-descent w.r.t. b in line
            // with thesis).
            //
            // temp = e_{tau FB} - s_{tau alpha}^T e_{tau FN}
            //

            #ifdef DEBUG
            std::cerr << "starting singular step\n";
            #endif

            start  = (problem.N_C)+1;
            finish = (problem.N)-nbad+1;

            #ifdef DO__FLOPS
            REG_n_step_sing;
            #endif

            d_alpha_pivot[(problem.N)-nbad+1-1] = 1.0;

            if ( nbad < problem.N_F )
            {
                //
                // Calculate the direction of zero curvature.
                //

                d_alpha_pivot.set_offset_start((problem.N_C)+1);
                d_alpha_pivot.set_offset_end((problem.N)-nbad);

                problem.opt_near_invert_factor(d_alpha_pivot,d_b);

                d_alpha_pivot.negate_it();
                d_b = -d_b;

                d_alpha_pivot.reset_offsets();
            }

            //
            // Work out the step scale and direction.
            //

            d_alpha_pivot.set_offset_start((problem.N_C)+1);
            d_alpha_pivot.set_offset_end((problem.N)-nbad+1);
            (problem.opt_get_e_free_pivot()).set_offset_start(1);
            (problem.opt_get_e_free_pivot()).set_offset_end((problem.N_F)-nbad+1);

            temp = ( d_alpha_pivot * (problem.opt_get_e_free_pivot()) );

            (problem.opt_get_e_free_pivot()).reset_offsets();
            d_alpha_pivot.reset_offsets();

            if ( problem.opt_get_tau_pivot((problem.N)-nbad) > 0 )
            {
                if ( temp >= 0.0 )
                {
                    theta = -1.1 * (problem.opt_get_h_pivot((problem.N)-nbad));
                }

                else
                {
                    theta = 1.1 * (problem.opt_get_h_pivot((problem.N)-nbad));
                }
            }

            else
            {
                if ( temp >= 0.0 )
                {
                    theta = -1.1 * -(problem.opt_get_v_pivot((problem.N)-nbad));
                }

                else
                {
                    theta = 1.1 * -(problem.opt_get_v_pivot((problem.N)-nbad));
                }
            }

            //
            // Scale and direct step.
            //

            d_alpha_pivot.set_offset_start((problem.N_C)+1);
            d_alpha_pivot.set_offset_end((problem.N)-nbad+1);

            d_alpha_pivot *= theta;
            d_b           *= theta;

            d_alpha_pivot.reset_offsets();

            if ( nbad >= 2 )
            {
                //
                // Fixup n_zero if necessary.
                //
        
                if ( n_zero > (problem.N_F)-nbad+1 )
                {
                    n_zero = (problem.N_F)-nbad+1;
                }
            }
        }
    }

    return;
}

void calc_step_scale_active(SVdata &problem, long &p, int &p_indic, fVECTOR &d_alpha_pivot, L_DOUBLE &d_b, fVECTOR &d_e_pivot, L_DOUBLE &d_f, long start, long finish)
{
    long i;

    L_DOUBLE ceta;
    L_DOUBLE deta;

    deta    = 1.0;
    p       = 0;
    p_indic = 0;

    if ( start <= finish )
    {
        //
        // Optimization note: use <= when checking a violation against the
        // worst violation so far to favour elements toward the end of V,
        // thus keeping e_zero elements and also speeding up cholesky
        // updates.
        //

        for ( i = start ; i <= finish ; i++ )
        {
            //
            // Note: A small step in the wrong direction can cause instant
            //       scaling of step to zero, which will result in an
            //       infinite loop.  The following statement is intended
            //       to prevent this occurence.
            //

// CAREFUL: if the tolerance is too large wrt G then troubles result.
//          Why?  Consider the case Gij very large.  Then there is a very
//          small change in alpha.  Ignoring the magnitude of Gij, it would
//          be tempting to say that the change in alpha is insignificant and
//          zero it.  However, Gij.d_alpha_i would be significant, so the
//          rounding of d_alpha_i will lead to a significant change in the
//          gradient (G.alpha component) and subsequently large errors.

            if ( ( d_alpha_pivot[i-1] > -ZZERO/10 ) && ( d_alpha_pivot[i-1] < ZZERO/10 ) )
            {
                d_alpha_pivot[i-1] = 0;
            }

            if ( d_alpha_pivot[i-1] < 0.0 )
            {
                //
                // Step is decreasing alpha[i-1]
                //

                if ( problem.opt_get_tau_pivot(i-1) > 0 )
                {
                    //
                    // The step is toward 0 from some positive value
                    //

                    ceta = -(problem.opt_get_alpha_pivot(i-1)) / d_alpha_pivot[i-1];

                    if ( ceta <= deta )
                    {
                        //
                        // Detected a constraint violation by this step,
                        // which is larger than any previous violation.
                        //

                        deta    =  ceta;
                        p       =  i;
                        p_indic =  1;
                    }
                }

                else
                {
                    //
                    // The step is toward v from some -ve value
                    //

                    ceta = ( (problem.opt_get_v_pivot(i-1)) - (problem.opt_get_alpha_pivot(i-1)) ) / d_alpha_pivot[i-1];

                    if ( ceta <= deta )
                    {
                        //
                        // Detected a constraint violation by this step,
                        // which is larger than any previous violation.
                        //

                        deta    =  ceta;
                        p       =  i;
                        p_indic = -1;
                    }
                }
            }

            else if ( d_alpha_pivot[i-1] > 0.0 )
            {
                //
                // Step is increasing alpha[i-1]
                //

                if ( problem.opt_get_tau_pivot(i-1) > 0 )
                {
                    //
                    // The step is toward h from some positive value
                    //
            
                    ceta = ( (problem.opt_get_h_pivot(i-1)) - (problem.opt_get_alpha_pivot(i-1)) ) / d_alpha_pivot[i-1];

                    if ( ceta <= deta )
                    {
                        //
                        // Detected a constraint violation by this step,
                        // which is larger than any previous violation.
                        //

                        deta    =  ceta;
                        p       =  i;
                        p_indic =  2;
                    }
                }

                else
                {
                    //
                    // The step is toward 0 from some negative value
                    //

                    ceta = -(problem.opt_get_alpha_pivot(i-1)) / d_alpha_pivot[i-1];

                    if ( ceta <= deta )
                    {
                        //
                        // Detected a constraint violation by this step,
                        // which is larger than any previous violation.
                        //

                        deta    =  ceta;
                        p       =  i;
                        p_indic = -2;
                    }
                }
            }
        }
    }

    if ( p_indic )
    {
        d_alpha_pivot.set_offset_start(start);
        d_alpha_pivot.set_offset_end(finish);
        d_e_pivot.set_offset_start(start);
        d_e_pivot.set_offset_end(finish);

        d_alpha_pivot *= deta;
        d_e_pivot     *= deta;
        d_b *= deta;
        d_f *= deta;

        d_e_pivot.reset_offsets();
        d_alpha_pivot.reset_offsets();
    }

    return;
}

void take_step_active(SVdata &problem, long &n_zero, int &f_zero, int p_indic, fVECTOR &d_alpha_pivot, L_DOUBLE &d_b, fVECTOR &d_e_pivot, L_DOUBLE &d_f, long start, long finish)
{
    if ( problem.N_F > 0 )
    {
        //
        // First, update f_zero
        //

        if ( p_indic == 0 )
        {
            f_zero = 1;
        }

        //
        // Update gradients for constrained variables
        //

        d_alpha_pivot.set_offset_start(start);
        d_alpha_pivot.set_offset_end(finish);
        d_e_pivot.set_offset_start(start);
        d_e_pivot.set_offset_end(finish);

        problem.opt_step_general_pivot(start,finish,d_alpha_pivot,d_b,d_e_pivot,d_f,n_zero);

        d_e_pivot.reset_offsets();
        d_alpha_pivot.reset_offsets();

        //
        // Finally, update n_zero
        //

        if ( p_indic == 0 )
        {
            n_zero = problem.N_F;
        }
    }

    return;
}

int test_optimality_active(SVdata &problem, long &q, int &q_indic, int f_zero, int p_indic, L_DOUBLE &sol_tol)
{
    long i;
    L_DOUBLE temp;
    L_DOUBLE lopt;

    if ( p_indic != 0 )
    {
        return 0;
    }

    //
    // NB: - lopt starts at -sol_tol to prevent non-termination due to
    //       numerical issues.
    //

    lopt    = -sol_tol;
    q       = 0;
    q_indic = 0;

    if ( problem.N_Z > 0 )
    {
        for ( i = 1 ; i <= problem.N_Z ; i++ )
        {
            //
            // Test to see if increasing alpha from 0 (if possible) is more
            // non-optimal than the least non-optimal found so far.
            //
            // Note: - the condition of least non-optimality may be dropped
            //         if f != 0 and this is the first possible change.  This
            //         basically means we can go up-hill if necessary, which
            //         prevents a possible infinite loop situation due to
            //         numerical problems.
            //

            temp = problem.opt_get_e_pos_pivot(i-1);

            #ifdef DEBUG
            std::cerr << "posgrad " << i << " = " << temp << "\n";
            #endif

            if ( ( ( temp <= lopt ) || ( ( f_zero == 0 ) && ( q_indic == 0 ) ) ) &&
                 ( ( (problem.opt_get_contype_pivot(i-1)) == 2 ) || ( (problem.opt_get_contype_pivot(i-1)) == 3 ) ) &&
                 ( ( (problem.N_F) > 0 ) || ( problem.f <= 0.0 ) || ( f_zero == 1 ) ) )
            {
                lopt    =  temp;
                q       =  i;
                q_indic =  1;
            }

            //
            // ... or should we decrease alpha below 0?
            //

            temp = -(problem.opt_get_e_neg_pivot(i-1));

            #ifdef DEBUG
            std::cerr << "neggrad " << i << " = " << temp << "\n";
            #endif

            if ( ( ( temp <= lopt ) || ( ( f_zero == 0 ) && ( q_indic == 0 ) ) ) &&
                 ( ( (problem.opt_get_contype_pivot(i-1)) == 1 ) || ( (problem.opt_get_contype_pivot(i-1)) == 3 ) ) &&
                 ( ( (problem.N_F) > 0 ) || ( problem.f >= 0.0 ) || ( f_zero == 1 ) ) )
            {
                lopt    =  temp;
                q       =  i;
                q_indic = -1;
            }
        }
    }

    if ( problem.N_L > 0 )
    {
        for ( i = (problem.N_Z)+1 ; i <= (problem.N_Z)+(problem.N_L) ; i++ )
        {
            //
            // Consider increasing alpha from v
            //

            temp =  problem.opt_get_e_neg_pivot(i-1);

            if ( ( ( temp <= lopt ) || ( ( f_zero == 0 ) && ( q_indic == 0 ) )      ) &&
                 ( ( (problem.N_F) > 0 ) || ( problem.f <= 0.0 ) || ( f_zero == 1 ) )    )
            {
                lopt    =  temp;
                q       =  i;
                q_indic = -2;
            }
        }
    }

    if ( problem.N_U > 0 )
    {
        for ( i = (problem.N_Z)+(problem.N_L)+1 ; i <= (problem.N_Z)+(problem.N_L)+(problem.N_U) ; i++ )
        {
            //
            // Consider decreasing alpha from h
            //

            temp = -(problem.opt_get_e_pos_pivot(i-1));

            if ( ( ( temp <= lopt ) || ( ( f_zero == 0 ) && ( q_indic == 0 ) )      ) &&
                 ( ( (problem.N_F) > 0 ) || ( problem.f >= 0.0 ) || ( f_zero == 1 ) )    )
            {
                lopt    =  temp;
                q       =  i;
                q_indic =  2;
            }
        }
    }

    if ( q_indic != 0 )
    {
        return 0;
    }

    return 1;
}

void fix_active_set_active(SVdata &problem, long p, int p_indic, long q, int q_indic, long &n_zero)
{
    //
    // First: if a constraint violation was detected when scaling
    //        the step, then constrain the relevant variable.
    //

    if ( p_indic )
    {
        if ( p-((problem.N_Z)+(problem.N_L)+(problem.N_U)) <= n_zero )
        {
            n_zero--;
        }

        switch ( p_indic )
        {
            case 1:
            {
                #ifdef DO__FLOPS
                REG_n_c_lower;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e1\n";
                std::cerr << "p = " << p << "\n";
                std::cerr << "dest = " << (problem.N_Z)+1 << "\n";
                #endif

                //
                // Constrain at zero from above
                //

                problem.opt_constrain_Z_pivot(p-1);

                break;
            }

            case 2:
            {
                #ifdef DO__FLOPS
                REG_n_c_upper;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e2\n";
                #endif

                //
                // Constrain at h from below
                //

                problem.opt_constrain_U_pivot(p-1);

                break;
            }

            case -2:
            {
                #ifdef DO__FLOPS
                REG_n_c_lower;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e3\n";
                #endif

                //
                // Constrain at 0 from below
                //

                problem.opt_constrain_Z_pivot(p-1);

                break;
            }

            case -1:
            {
                #ifdef DO__FLOPS
                REG_n_c_upper;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e4\n";
                #endif

                //
                // Constrain at v from above
                //

                problem.opt_constrain_L_pivot(p-1);

                break;
            }

            default:
            {
                L_THROW(0);

                break;
            }
        }
    }

    else
    {
        //
        // Second: Unconstrain largest inoptimality found.
        //

        switch ( q_indic )
        {
            case 1:
            {
                #ifdef DO__FLOPS
                REG_n_u_lower;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e5\n";
                #endif

                //
                // unconstrain zero to above
                //

                problem.opt_free_U_pivot(q-1);

                break;
            }

            case 2:
            {
                #ifdef DO__FLOPS
                REG_n_u_upper;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e6\n";
                #endif

                //
                // unconstrain upper bound h
                //

                problem.opt_free_U_pivot(q-1);

                break;
            }

            case -1:
            {
                #ifdef DO__FLOPS
                REG_n_u_lower;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e7\n";
                #endif

                //
                // unconstrain zero to below
                //

                problem.opt_free_L_pivot(q-1);

                break;
            }

            case -2:
            {
                #ifdef DO__FLOPS
                REG_n_u_upper;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e8\n";
                #endif

                //
                // unconstrain lower bound v
                //

                problem.opt_free_L_pivot(q-1);

                break;
            }

            default:
            {
                L_THROW(0);

                break;
            }
        }
    }

    return;

    q = 0;
}


unsigned long solve_active_fixed_bias(SVdata &problem, L_DOUBLE &sol_tol, unsigned long epochs, volatile int *async_exit_flag, volatile unsigned long *loop_count)
{
    long i;

    long start;
    long finish;

    long n_zero;

    fVECTOR d_alpha_pivot('r',problem.N);
    fVECTOR d_e_pivot('r',problem.N);

    long p = 0;
    int p_indic = 0;

    long q = 0;
    int q_indic = 0;

    //
    // Set initial state for n_zero
    //

    n_zero = 0;

    if ( problem.N_F > 0 )
    {
        for ( i = (problem.N_C)+1 ; i <= problem.N ; i++ )
        {
            if ( problem.opt_get_e_pivot(i-1) != 0.0 )
            {
                break;
            }

            n_zero++;
        }
    }

    //
    // Main algorithmic loop.
    //

    for ( ; ( ( (*loop_count) <= epochs ) && ( !(*async_exit_flag) ) ) ; (*loop_count)++ )
    {
        #ifdef DEBUG
        std::cerr << "===========" << (*loop_count) << "=============\n";
        std::cerr << problem << "\n\n";
        #endif

        #ifdef DO__FLOPS
        REG_n_iterations;
        #endif

        //
        // Calculate the step (ignoring bounds on free variables).  This
        // will also detect (and fix) any changes to n_zero.
        //

        start  = (problem.N_C)+1;
        finish = problem.N;

        #ifdef DEBUG
        std::cerr << "phantom step a\n";
        #endif

        calc_step_active_fixed_bias(problem,n_zero,d_alpha_pivot,d_e_pivot,start,finish);

        //
        // Calculate scale factor for step.  The scale factor is deta, and
        // p and p_indic mark the alpha which caused this as described
        // previously.
        //

        #ifdef DEBUG
        std::cerr << "Phantom step b\n";
        #endif

        calc_step_scale_active_fixed_bias(problem,p,p_indic,d_alpha_pivot,d_e_pivot,start,finish);

        //
        // Take scaled step, fix f_zero and n_zero appropriately.
        //

        #ifdef DEBUG
        std::cerr << "Phantom step c\n";
        std::cerr << "d_alpha: " << d_alpha_pivot << "\n\n";
        std::cerr << "d_e: " << d_e_pivot << "\n\n";
        #endif

        take_step_active_fixed_bias(problem,n_zero,p_indic,d_alpha_pivot,d_e_pivot,start,finish);

        #ifdef DEBUG
        std::cerr << "===========" << (*loop_count) << " phantom=============\n";
        std::cerr << "d_alpha: " << d_alpha_pivot << "\n\n";
        std::cerr << "d_e: " << d_e_pivot << "\n\n";
        std::cerr << "start = " << start << "\n";
        std::cerr << "finish = " << finish << "\n";
        std::cerr << problem << "\n\n";
        #endif

        //
        // Check optimality of solution, exit if optimal.  Also find the
        // most inoptimal variable if solution is not optimal.
        //

        #ifdef DEBUG
        std::cerr << "Phantom step d\n";
        #endif

        if ( test_optimality_active_fixed_bias(problem,q,q_indic,p_indic,sol_tol) )
        {
            break;
        }

        //
        // Modify the active set.  Fix n_zero if required.
        //

        #ifdef DEBUG
        std::cerr << "Phantom step e\n";
        #endif

        fix_active_set_active_fixed_bias(problem,p,p_indic,q,q_indic,n_zero);
    }

    return (*loop_count);
}

void calc_step_active_fixed_bias(SVdata &problem, long &n_zero, fVECTOR &d_alpha_pivot, fVECTOR &d_e_pivot, long &start, long &finish)
{
    long nbad;

    L_DOUBLE theta;
    L_DOUBLE temp;

    //
    // Default step is zero
    //

    d_e_pivot     = 0.0;
    d_alpha_pivot = 0.0;

    //
    // We only calculate a step if N_F > 0
    //

    if ( problem.N_F > 0 )
    {
        //
        // Find how much (if any) of V is singular
        //
        // Note that if fixed bias is selected then V functions will ensure
        // that d_b is set appropriately to zero.
        //

        if ( ( nbad = problem.opt_get_nbad_factor() ) == 0 )
        {
            //
            // V is non-singular, so step is Newton.
            //

            start  = (problem.N_C)+1;
            finish = problem.N;

            #ifdef DO__FLOPS
            REG_n_step_nonsing;
            #endif

            //
            // If e is all zero, the step is trivial so no need to bother
            // actually calculating it.
            //

            if ( n_zero < problem.N_F )
            {
                d_e_pivot.set_offset_start((problem.N_C)+1);
                d_e_pivot.set_offset_end(problem.N);
                d_alpha_pivot.set_offset_start((problem.N_C)+1);
                d_alpha_pivot.set_offset_end(problem.N);

                problem.opt_minverse_factor(d_alpha_pivot,n_zero,0);

                d_e_pivot = problem.opt_get_e_free_pivot();

                d_alpha_pivot.negate_it();
                d_e_pivot.negate_it();

                d_alpha_pivot.reset_offsets();
                d_e_pivot.reset_offsets();

                #ifdef DEBUG
                std::cerr << "phantom phantom\n";
                std::cerr << start << " - " << finish << "\n";
                std::cerr << (problem.N_C) << " - " << problem.N << "\n";
                std::cerr << "d_e: " << d_e_pivot << "\n\n";
                std::cerr << "d_alpha: " << d_alpha_pivot << "\n\n";
                #endif
            }
        }
    
        else
        {
            //
            // V is singular, so we step in a direction of linear non-ascent
            // w.r.t. alpha (ignoring possible non-descent w.r.t. b in line
            // with thesis).
            //
            // temp = e_{tau FB} - s_{tau alpha}^T e_{tau FN}
            //

            start  = (problem.N_C)+1;
            finish = (problem.N)-nbad+1;

            #ifdef DO__FLOPS
            REG_n_step_sing;
            #endif

            d_alpha_pivot[(problem.N)-nbad+1-1] = 1.0;

            if ( nbad < problem.N_F )
            {
                //
                // Calculate the direction of zero curvature.
                //

                d_alpha_pivot.set_offset_start((problem.N_C)+1);
                d_alpha_pivot.set_offset_end((problem.N)-nbad);

                problem.opt_near_invert_factor(d_alpha_pivot);

                d_alpha_pivot.negate_it();

                d_alpha_pivot.reset_offsets();
            }

            //
            // Work out the step scale and direction.
            //

            d_alpha_pivot.set_offset_start((problem.N_C)+1);
            d_alpha_pivot.set_offset_end((problem.N)-nbad+1);
            (problem.opt_get_e_free_pivot()).set_offset_start(1);
            (problem.opt_get_e_free_pivot()).set_offset_end((problem.N_F)-nbad+1);

            temp = ( d_alpha_pivot * (problem.opt_get_e_free_pivot()) );

            (problem.opt_get_e_free_pivot()).reset_offsets();
            d_alpha_pivot.reset_offsets();

            if ( problem.opt_get_tau_pivot((problem.N)-nbad) > 0 )
            {
                if ( temp >= 0.0 )
                {
                    theta = -1.1 * (problem.opt_get_h_pivot((problem.N)-nbad));
                }

                else
                {
                    theta = 1.1 * (problem.opt_get_h_pivot((problem.N)-nbad));
                }
            }

            else
            {
                if ( temp >= 0.0 )
                {
                    theta = -1.1 * -(problem.opt_get_v_pivot((problem.N)-nbad));
                }

                else
                {
                    theta = 1.1 * -(problem.opt_get_v_pivot((problem.N)-nbad));
                }
            }

            //
            // Scale and direct step.
            //

            d_alpha_pivot.set_offset_start((problem.N_C)+1);
            d_alpha_pivot.set_offset_end((problem.N)-nbad+1);

            d_alpha_pivot *= theta;

            d_alpha_pivot.reset_offsets();

            if ( nbad >= 2 )
            {
                //
                // Fixup n_zero if necessary.
                //
        
                if ( n_zero > (problem.N_F)-nbad+1 )
                {
                    n_zero = (problem.N_F)-nbad+1;
                }
            }
        }
    }

    return;
}

void calc_step_scale_active_fixed_bias(SVdata &problem, long &p, int &p_indic, fVECTOR &d_alpha_pivot, fVECTOR &d_e_pivot, long start, long finish)
{
    long i;

    L_DOUBLE ceta;
    L_DOUBLE deta;

    deta    = 1.0;
    p       = 0;
    p_indic = 0;

    if ( start <= finish )
    {
        //
        // Optimization note: use <= when checking a violation against the
        // worst violation so far to favour elements toward the end of V,
        // thus keeping e_zero elements and also speeding up cholesky
        // updates.
        //

        for ( i = start ; i <= finish ; i++ )
        {
            //
            // Note: A small step in the wrong direction can cause instant
            //       scaling of step to zero, which will result in an
            //       infinite loop.  The following statement is intended
            //       to prevent this occurence.
            //

            if ( ( d_alpha_pivot[i-1] > -ZZERO/10 ) && ( d_alpha_pivot[i-1] < ZZERO/10 ) )
            {
                d_alpha_pivot[i-1] = 0;
            }

            if ( d_alpha_pivot[i-1] < 0.0 )
            {
                //
                // Step is decreasing alpha[i-1]
                //

                if ( problem.opt_get_tau_pivot(i-1) > 0 )
                {
                    //
                    // The step is toward 0 from some positive value
                    //

                    ceta = -(problem.opt_get_alpha_pivot(i-1)) / d_alpha_pivot[i-1];

                    if ( ceta <= deta )
                    {
                        //
                        // Detected a constraint violation by this step,
                        // which is larger than any previous violation.
                        //

                        deta    =  ceta;
                        p       =  i;
                        p_indic =  1;
                    }
                }

                else
                {
                    //
                    // The step is toward v from some -ve value
                    //

                    ceta = ( (problem.opt_get_v_pivot(i-1)) - (problem.opt_get_alpha_pivot(i-1)) ) / d_alpha_pivot[i-1];

                    if ( ceta <= deta )
                    {
                        //
                        // Detected a constraint violation by this step,
                        // which is larger than any previous violation.
                        //

                        deta    =  ceta;
                        p       =  i;
                        p_indic = -1;
                    }
                }
            }

            else if ( d_alpha_pivot[i-1] > 0.0 )
            {
                //
                // Step is increasing alpha[i-1]
                //

                if ( problem.opt_get_tau_pivot(i-1) > 0 )
                {
                    //
                    // The step is toward h from some positive value
                    //
            
                    ceta = ( (problem.opt_get_h_pivot(i-1)) - (problem.opt_get_alpha_pivot(i-1)) ) / d_alpha_pivot[i-1];

                    if ( ceta <= deta )
                    {
                        //
                        // Detected a constraint violation by this step,
                        // which is larger than any previous violation.
                        //

                        deta    =  ceta;
                        p       =  i;
                        p_indic =  2;
                    }
                }

                else
                {
                    //
                    // The step is toward 0 from some negative value
                    //

                    ceta = -(problem.opt_get_alpha_pivot(i-1)) / d_alpha_pivot[i-1];

                    if ( ceta <= deta )
                    {
                        //
                        // Detected a constraint violation by this step,
                        // which is larger than any previous violation.
                        //

                        deta    =  ceta;
                        p       =  i;
                        p_indic = -2;
                    }
                }
            }
        }
    }

    if ( p_indic )
    {
        d_alpha_pivot.set_offset_start(start);
        d_alpha_pivot.set_offset_end(finish);
        d_e_pivot.set_offset_start(start);
        d_e_pivot.set_offset_end(finish);

        d_alpha_pivot *= deta;
        d_e_pivot     *= deta;

        d_e_pivot.reset_offsets();
        d_alpha_pivot.reset_offsets();
    }

    return;
}

void take_step_active_fixed_bias(SVdata &problem, long &n_zero, int p_indic, fVECTOR &d_alpha_pivot, fVECTOR &d_e_pivot, long start, long finish)
{
    if ( problem.N_F > 0 )
    {
        //
        // Update gradients for constrained variables
        //

        d_alpha_pivot.set_offset_start(start);
        d_alpha_pivot.set_offset_end(finish);
        d_e_pivot.set_offset_start(start);
        d_e_pivot.set_offset_end(finish);

        problem.opt_step_general_pivot(start,finish,d_alpha_pivot,0.0,d_e_pivot,0.0,n_zero);

        d_e_pivot.reset_offsets();
        d_alpha_pivot.reset_offsets();

        //
        // Update n_zero.
        //

        if ( p_indic == 0 )
        {
            n_zero = problem.N_F;
        }
    }

    return;
}

int test_optimality_active_fixed_bias(SVdata &problem, long &q, int &q_indic, int p_indic, L_DOUBLE &sol_tol)
{
    long i;
    L_DOUBLE temp;
    L_DOUBLE lopt;

    if ( p_indic != 0 )
    {
        return 0;
    }

    //
    // NB: - lopt starts at -sol_tol to prevent non-termination due to
    //       numerical issues.
    //

    lopt    = -sol_tol;
    q       = 0;
    q_indic = 0;

    if ( problem.N_Z > 0 )
    {
        for ( i = 1 ; i <= problem.N_Z ; i++ )
        {
            //
            // Test to see if increasing alpha from 0 (if possible) is more
            // non-optimal than the least non-optimal found so far.
            //

            temp = problem.opt_get_e_pos_pivot(i-1);

            if ( ( temp <= lopt ) && ( ( (problem.opt_get_contype_pivot(i-1)) == 2 ) || ( (problem.opt_get_contype_pivot(i-1)) == 3 ) ) )
            {
                lopt    =  temp;
                q       =  i;
                q_indic =  1;
            }

            //
            // ... or should we decrease alpha below 0?
            //

            temp = -(problem.opt_get_e_neg_pivot(i-1));

            if ( ( temp <= lopt ) && ( ( (problem.opt_get_contype_pivot(i-1)) == 1 ) || ( (problem.opt_get_contype_pivot(i-1)) == 3 ) ) )
            {
                lopt    =  temp;
                q       =  i;
                q_indic = -1;
            }
        }
    }

    if ( problem.N_L > 0 )
    {
        for ( i = (problem.N_Z)+1 ; i <= (problem.N_Z)+(problem.N_L) ; i++ )
        {
            //
            // Consider increasing alpha from v
            //

            temp =  problem.opt_get_e_neg_pivot(i-1);

            if ( temp <= lopt )
            {
                lopt    =  temp;
                q       =  i;
                q_indic = -2;
            }
        }
    }

    if ( problem.N_U > 0 )
    {
        for ( i = (problem.N_Z)+(problem.N_L)+1 ; i <= (problem.N_Z)+(problem.N_L)+(problem.N_U) ; i++ )
        {
            //
            // Consider decreasing alpha from h
            //

            temp = -(problem.opt_get_e_pos_pivot(i-1));

            if ( temp <= lopt )
            {
                lopt    =  temp;
                q       =  i;
                q_indic =  2;
            }
        }
    }

    if ( q_indic != 0 )
    {
        return 0;
    }

    return 1;
}

void fix_active_set_active_fixed_bias(SVdata &problem, long p, int p_indic, long q, int q_indic, long &n_zero)
{
    //
    // First: if a constraint violation was detected when scaling
    //        the step, then constrain the relevant variable.
    //

    if ( p_indic )
    {
        if ( p-((problem.N_Z)+(problem.N_L)+(problem.N_U)) <= n_zero )
        {
            n_zero--;
        }

        switch ( p_indic )
        {
            case 1:
            {
                #ifdef DO__FLOPS
                REG_n_c_lower;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e1\n";
                std::cerr << "p = " << p << "\n";
                std::cerr << "dest = " << (problem.N_Z)+1 << "\n";
                #endif

                //
                // Constrain at zero from above
                //

                problem.opt_constrain_Z_pivot(p-1);

                break;
            }

            case 2:
            {
                #ifdef DO__FLOPS
                REG_n_c_upper;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e2\n";
                #endif

                //
                // Constrain at h from below
                //

                problem.opt_constrain_U_pivot(p-1);

                break;
            }

            case -2:
            {
                #ifdef DO__FLOPS
                REG_n_c_lower;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e3\n";
                #endif

                //
                // Constrain at 0 from below
                //

                problem.opt_constrain_Z_pivot(p-1);

                break;
            }

            case -1:
            {
                #ifdef DO__FLOPS
                REG_n_c_upper;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e4\n";
                #endif

                //
                // Constrain at v from above
                //

                problem.opt_constrain_L_pivot(p-1);

                break;
            }

            default:
            {
                L_THROW(0);

                break;
            }
        }
    }

    else
    {
        //
        // Second: Unconstrain largest inoptimality found.
        //

        switch ( q_indic )
        {
            case 1:
            {
                #ifdef DO__FLOPS
                REG_n_u_lower;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e5\n";
                #endif

                //
                // unconstrain zero to above
                //

                problem.opt_free_U_pivot(q-1);

                break;
            }

            case 2:
            {
                #ifdef DO__FLOPS
                REG_n_u_upper;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e6\n";
                #endif

                //
                // unconstrain upper bound h
                //

                problem.opt_free_U_pivot(q-1);

                break;
            }

            case -1:
            {
                #ifdef DO__FLOPS
                REG_n_u_lower;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e7\n";
                #endif

                //
                // unconstrain zero to below
                //

                problem.opt_free_L_pivot(q-1);

                break;
            }

            case -2:
            {
                #ifdef DO__FLOPS
                REG_n_u_upper;
                #endif

                #ifdef DEBUG
                std::cerr << "Phantom step e8\n";
                #endif

                //
                // unconstrain lower bound v
                //

                problem.opt_free_L_pivot(q-1);

                break;
            }

            default:
            {
                L_THROW(0);

                break;
            }
        }
    }

    return;

    q = 0;
}




/***********************************************************************

                        Platt SMO optimiser
                        ===================

***********************************************************************/


unsigned long solve_SMO(SVdata &problem, L_DOUBLE &sol_tol, unsigned long epochs, volatile int *async_exit_flag, volatile unsigned long *loop_count, JTYPE *J)
{
    int f_zero;
    long i;
    long numChanged = 0;
    int examineAll = 1;

    //
    // f_zero = 0 - if f != 0
    //        = 1 - if f == 0
    //

    if ( problem.f == 0.0 ) { f_zero = 1; }
    else                    { f_zero = 0; }

    //
    // Main algorithmic loop.
    //

    //
    // This is not quite Platt's original implementation, but it is
    // essentially the same.  In the singular case, we move in a
    // direction of linear descent rather then giving up (it can be
    // shown that such a direction will always exist if the error is
    // appropriate and f = 0).
    //
    
    #ifdef DEBUG
    std::cerr << "Begin SMO algorithm\n";
    std::cerr << problem << "\n";
    std::cerr << "=====================\n\n\n";
    #endif

    for ( ; ( ( (*loop_count) <= epochs ) && ( !(*async_exit_flag) ) ) ; )
    {
        numChanged = 0;

        for ( i = 1 ; i <= problem.N ; i++ )
        {
            switch ( (problem.fast_tau)[i-1] )
            {
                case 0:
                {
                    /*
                       point constrained at zero.  If examineAll is set then
                       proceed to examine all posible alternative "states"
                       for this variable.
                    */

                    #ifdef DEBUG
                    std::cerr << "Examining 0 point\n";
                    #endif

                    if ( examineAll )
                    {
                        switch ( (problem.fast_contype)[i-1] )
                        {
                            case 0:
                            {
                                break;
                            }

                            case 1:
                            {
                                #ifdef DEBUG
                                std::cerr << "Trying free low\n";
                                #endif

                                #ifdef DEBUG
                                std::cerr << "...examining...\n";
                                #endif

                                numChanged += examineExample_SMO(i,-1,f_zero,problem,sol_tol,J);

                                break;
                            }

                            case 2:
                            {
                                #ifdef DEBUG
                                std::cerr << "Trying free upp\n";
                                #endif

                                #ifdef DEBUG
                                std::cerr << "...examining...\n";
                                #endif

                                numChanged += examineExample_SMO(i,+1,f_zero,problem,sol_tol,J);

                                break;
                            }

                            case 3:
                            {
                                if ( rand() % 2 )
                                {
                                    #ifdef DEBUG
                                    std::cerr << "Trying free low\n";
                                    #endif

                                    #ifdef DEBUG
                                    std::cerr << "...examining...\n";
                                    #endif

                                    if ( examineExample_SMO(i,-1,f_zero,problem,sol_tol,J) )
                                    {
                                        numChanged++;
                                    }

                                    else
                                    {
                                        #ifdef DEBUG
                                        std::cerr << "Trying free upp\n";
                                        #endif

                                        #ifdef DEBUG
                                        std::cerr << "...examining...\n";
                                        #endif

                                        numChanged += examineExample_SMO(i,+1,f_zero,problem,sol_tol,J);
                                    }
                                }

                                else
                                {
                                    #ifdef DEBUG
                                    std::cerr << "Trying free upp\n";
                                    #endif

                                    #ifdef DEBUG
                                    std::cerr << "...examining...\n";
                                    #endif

                                    if ( examineExample_SMO(i,+1,f_zero,problem,sol_tol,J) )
                                    {
                                        numChanged++;
                                    }

                                    else
                                    {
                                        #ifdef DEBUG
                                        std::cerr << "Trying free low\n";
                                        #endif

                                        #ifdef DEBUG
                                        std::cerr << "...examining...\n";
                                        #endif

                                        numChanged += examineExample_SMO(i,-1,f_zero,problem,sol_tol,J);
                                    }
                                }

                                break;
                            }

                            default:
                            {
                                L_THROW(0);

                                break;
                            }
                        }
                    }

                    break;
                }

                case -1:
                {
                    /*
                       Variable is free, so examine step
                    */

                    #ifdef DEBUG
                    std::cerr << "Examining +-1 point\n";
                    #endif

                    numChanged += examineExample_SMO(i,-1,f_zero,problem,sol_tol,J);

                    break;
                }

                case +1:
                {
                    /*
                       Variable is free, so examine step
                    */

                    #ifdef DEBUG
                    std::cerr << "Examining +-1 point\n";
                    #endif

                    numChanged += examineExample_SMO(i,+1,f_zero,problem,sol_tol,J);

                    break;
                }

                case +2:
                {
                    /*
                       point constrained at h.  If examineAll is set then
                       proceed to examine freeing this variable.
                    */

                    #ifdef DEBUG
                    std::cerr << "Examining +2 point\n";
                    #endif

                    if ( examineAll )
                    {
                        numChanged += examineExample_SMO(i,+1,f_zero,problem,sol_tol,J);
                    }

                    break;
                }

                case -2:
                {
                    /*
                       point constrained at v.  If examineAll is set then
                       proceed to examine freeing this variable.
                    */

                    #ifdef DEBUG
                    std::cerr << "Examining -2 point\n";
                    #endif

                    if ( examineAll )
                    {
                        numChanged += examineExample_SMO(i,-1,f_zero,problem,sol_tol,J);
                    }

                    break;
                }

                default:
                {
                    L_THROW(0);

                    break;
                }
            }
        }

        (*loop_count) += numChanged;

        if ( examineAll )
        {
            if ( numChanged == 0 )
            {
                goto getout;
            }

            examineAll = 0;
        }

        else if ( numChanged == 0 )
        {
            examineAll = 1;
        }
    }

    getout:

    return (*loop_count);
}

int examineExample_SMO(long i2, int tau2, int f_zero, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J)
{
    /*
       NB: this code has alpha <= 0 for d == -1 and alpha >= 0 for d == +1.
           However, SMO (and most other papers) use alpha >= 0 is both
           cases.
    */

    long i1start;
    long i1;
    L_DOUBLE E1;

    L_DOUBLE alph2;
    L_DOUBLE E2;
    L_DOUBLE r2;
    L_DOUBLE C2;

    int firstrunflag;

    alph2 = ( tau2 == +1 ) ? (problem.fast_alpha)[i2-1]  : -(problem.fast_alpha)[i2-1];
    E2    = ( tau2 == +1 ) ? problem.opt_get_e_pos(i2-1) : problem.opt_get_e_neg(i2-1);
    r2    = ( tau2 == +1 ) ? E2                          : -E2;
    C2    = ( tau2 == +1 ) ? (problem.fast_h)[i2-1]      : -(problem.fast_v)[i2-1];

    #ifdef DEBUG
    std::cerr << "\nExamine Example\n===============\n";
    std::cerr << "i2    = " << i2 << "\n";
    std::cerr << "tau2  = " << tau2 << "\n";
    std::cerr << "alph2 = " << alph2 << "\n";
    std::cerr << "E2    = " << E2 << "\n";
    std::cerr << "r2    = " << r2 << "\n";
    std::cerr << "C2    = " << C2 << "\n";
    #endif

    if ( ( ( r2 < -sol_tol ) && ( alph2 < C2 ) ) || ( ( r2 > sol_tol ) && ( alph2 > 0.0 ) ) )
    {
        #ifdef DEBUG
        std::cerr << "Second choice heuristic\n";
        #endif

        if ( ( i1 = secondChoice_SMO(E1,E2,problem) ) < 0 )
        {
            if ( takeStep_SMO(-i1,i2,-1,tau2,E1,E2,f_zero,problem,sol_tol,J) )
            {
                return 1;
            }
        }

        else if ( i1 > 0 )
        {
            if ( takeStep_SMO(i1,i2,+1,tau2,E1,E2,f_zero,problem,sol_tol,J) )
            {
                return 1;
            }

        }

        #ifdef DEBUG
        std::cerr << "Second choice first fallback\n";
        #endif

        i1start = ( rand() % (problem.N) ) + 1;
        firstrunflag = 1;

        for ( i1 = i1start ; ( i1 != i1start ) || firstrunflag ; ( i1 = ( ( i1 == (problem.N) ) ? 1 : (i1+1) ) ) )
        {
            #ifdef DEBUG
            std::cerr << "...trying" << i1 << "...\n";
            #endif

            if ( (problem.fast_tau)[i1-1] == +1 )
            {
                E1 = problem.opt_get_e_pos(i1-1);

                if ( takeStep_SMO(i1,i2,+1,tau2,E1,E2,f_zero,problem,sol_tol,J) )
                {
                    return 1;
                }
            }

            else if ( (problem.fast_tau)[i1-1] == -1 )
            {
                E1 = problem.opt_get_e_neg(i1-1);

                if ( takeStep_SMO(i1,i2,-1,tau2,E1,E2,f_zero,problem,sol_tol,J) )
                {
                    return 1;
                }
            }

            firstrunflag = 0;
        }

        #ifdef DEBUG
        std::cerr << "Second choise second fallback\n";
        #endif

        i1start = ( rand() % (problem.N) ) + 1;
        firstrunflag = 1;

        for ( i1 = i1start ; ( i1 != i1start ) || firstrunflag ; ( i1 = ( ( i1 == (problem.N) ) ? 1 : (i1+1) ) ) )
        {
            #ifdef DEBUG
            std::cerr << "...trying" << i1 << "...\n";
            #endif

            switch ( (problem.fast_tau)[i1-1] )
            {
                case 0:
                {
                    switch ( (problem.fast_contype)[i1-1] )
                    {
                        case 0:
                        {
                            break;
                        }

                        case 1:
                        {
                            #ifdef DEBUG
                            std::cerr << "Second choice free low " << i1 << "\n";
                            #endif

                            E1 = problem.opt_get_e_neg(i1-1);

                            if ( takeStep_SMO(i1,i2,-1,tau2,E1,E2,f_zero,problem,sol_tol,J) )
                            {
                                return 1;
                            }

                            break;
                        }

                        case 2:
                        {
                            #ifdef DEBUG
                            std::cerr << "Second choice free upp " << i1 << "\n";
                            #endif

                            E1 = problem.opt_get_e_pos(i1-1);

                            if ( takeStep_SMO(i1,i2,+1,tau2,E1,E2,f_zero,problem,sol_tol,J) )
                            {
                                return 1;
                            }

                            break;
                        }

                        case 3:
                        {
                            if ( rand() % 2 )
                            {
                                #ifdef DEBUG
                                std::cerr << "Second choice free low " << i1 << "\n";
                                #endif

                                E1 = problem.opt_get_e_neg(i1-1);

                                if ( takeStep_SMO(i1,i2,-1,tau2,E1,E2,f_zero,problem,sol_tol,J) )
                                {
                                    return 1;
                                }

                                else
                                {
                                    #ifdef DEBUG
                                    std::cerr << "Second choice free upp " << i1 << "\n";
                                    #endif

                                    E1 = problem.opt_get_e_pos(i1-1);

                                    if ( takeStep_SMO(i1,i2,+1,tau2,E1,E2,f_zero,problem,sol_tol,J) )
                                    {
                                        return 1;
                                    }
                                }
                            }

                            else
                            {
                                #ifdef DEBUG
                                std::cerr << "Second choice free upp " << i1 << "\n";
                                #endif

                                E1 = problem.opt_get_e_pos(i1-1);

                                if ( takeStep_SMO(i1,i2,+1,tau2,E1,E2,f_zero,problem,sol_tol,J) )
                                {
                                    return 1;
                                }

                                else
                                {
                                    #ifdef DEBUG
                                    std::cerr << "Second choice free low " << i1 << "\n";
                                    #endif

                                    E1 = problem.opt_get_e_neg(i1-1);

                                    if ( takeStep_SMO(i1,i2,-1,tau2,E1,E2,f_zero,problem,sol_tol,J) )
                                    {
                                        return 1;
                                    }
                                }
                            }

                            break;
                        }

                        default:
                        {
                            L_THROW(0);

                            break;
                        }
                    }

                    break;
                }

                case +1:
                case +2:
                {
                    E1 = problem.opt_get_e_pos(i1-1);

                    if ( takeStep_SMO(i1,i2,+1,tau2,E1,E2,f_zero,problem,sol_tol,J) )
                    {
                        return 1;
                    }

                    break;
                }

                case -1:
                case -2:
                {
                    E1 = problem.opt_get_e_neg(i1-1);

                    if ( takeStep_SMO(i1,i2,-1,tau2,E1,E2,f_zero,problem,sol_tol,J) )
                    {
                        return 1;
                    }

                    break;
                }

                default:
                {
                    L_THROW(0);

                    break;
                }
            }

            firstrunflag = 0;
        }
    }

    return 0;
}

long secondChoice_SMO(L_DOUBLE &E1, L_DOUBLE &E2, SVdata &problem)
{
    long istart;
    long i;
    int firstTry = 1;
    int tau;
    L_DOUBLE E;

    long i1 = 0;
    int tau1 = 0;

    E1 = 0.0;

    istart = ( rand() % (problem.N) ) + 1;

    for ( i = ( ( istart == (problem.N) ) ? 1 : (istart+1) ) ; i != istart ; ( i = ( ( i == (problem.N) ) ? 1 : (i+1) ) ) )
    {
        tau = (problem.fast_tau)[i-1];

        if ( tau == +1 )
        {
            E = problem.opt_get_e_pos(i-1);

            if ( firstTry || ( ( E2 >= 0.0 ) && ( E < E1 ) ) || ( ( E2 <= 0.0 ) && ( E > E1 ) ) )
            {
                i1   = i;
                E1   = E;
                tau1 = +1;
            }

            firstTry = 0;
        }

        else if ( tau == -1 )
        {
            E = problem.opt_get_e_neg(i-1);

            if ( firstTry || ( ( E2 >= 0.0 ) && ( E < E1 ) ) || ( ( E2 <= 0.0 ) && ( E > E1 ) ) )
            {
                i1   = i;
                E1   = E;
                tau1 = -1;
            }

            firstTry = 0;
        }
    }

    return tau1*i1;
}

int takeStep_SMO(long i1, long i2, int tau1, int tau2, L_DOUBLE &e1, L_DOUBLE &e2, int f_zero, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J)
{
    if ( i1 == i2 )
    {
        return 0;
    }

    L_DOUBLE d_alpha1 = 0.0;
    L_DOUBLE d_alpha2 = 0.0;
    L_DOUBLE d_b = 0.0;
    L_DOUBLE d_f = 0.0;
    L_DOUBLE d_J_epart = 0.0;
    int f_zero_next = f_zero;

    if ( trial_step_SMO(d_J_epart,i1,i2,tau1,tau2,e1,e2,d_alpha1,d_alpha2,d_b,d_f,f_zero,f_zero_next,problem,sol_tol) == 0 )
    {
        return 0;
    }

    return actually_take_step_SMO(d_J_epart,i1,i2,tau1,tau2,d_alpha1,d_alpha2,d_b,d_f,f_zero,f_zero_next,problem,sol_tol,J);
}

int trial_step_SMO(L_DOUBLE &d_J_epart, long i1, long i2, int tau1, int tau2, L_DOUBLE &e1, L_DOUBLE &e2, L_DOUBLE &d_alpha1, L_DOUBLE &d_alpha2, L_DOUBLE &d_b, L_DOUBLE &d_f, int &f_zero, int &f_zero_next, SVdata &problem, L_DOUBLE &sol_tol)
{
    #ifdef DEBUG
    std::cerr << "Try step " << i1 << " - " << i2 << "\n";
    #endif

    //
    // It may be observed that the SMO approach is simply a special case of
    // active set, where only 2 variables are free.  To this end, it will
    // be noted that:
    //
    // [ d_b      ]          [ 0  1   1  ]   [  f  ]
    // [ d_alpha1 ] = - inv( [ 1 K11 K12 ] ) [ e_1 ]
    // [ d_alpha2 ]          [ 1 K12 K22 ]   [ e_2 ]
    //
    //                       1            [  K11.K22 - K12.K12 ]
    //              = ----------------- ( [      K12 - K22     ] f
    //                K11 + K22 - 2.K12   [      K12 - K11     ]
    //
    //                  [ K12 - K22 ]       [ K12 - K11 ]
    //                + [    -1     ] e_1 + [    +1     ] e_2 )
    //                  [    +1     ]       [    -1     ]
    //
    // [ d_f  ]   [ -d_f  ]
    // [ d_e1 ] = [ -d_e1 ]
    // [ d_e2 ]   [ -d_e2 ]
    //
    // Of course, this is only true if the 2 hessian is non-singular - ie.
    // K11 + K22 - 2.K12 != 0.  If this does not hold, no problem - just
    // find a direction of linear non-ascent wrt alpha_1 and alpha_2,
    // ignoring f in this case.
    //
    // NB: we are ignoring the SMO method here, as it involves evaluating
    //     the entire goddam objective, which is tres-costly one would
    //     think!
    //
    // Specifically, in the singular case:
    //
    // [ d_b      ] = - theta inv( [ 0  1  ] ) [  1  ]
    // [ d_alpha1 ]                [ 1 K11 ]   [ K12 ]
    //
    //              = - theta (1/-1) [ K11 -1 ] [  1  ]
    //                               [ -1  0  ] [ K12 ]
    //
    //              = theta [ K11 -1 ] [  1  ]
    //                      [ -1  0  ] [ K12 ]
    //
    //              = theta [ K11 - K12 ]
    //                      [    -1     ]
    //
    // d_alpha2 = theta
    //
    // [ d_f  ]   [ 0 ]
    // [ d_e1 ] = [ 0 ]
    // [ d_e2 ]   [ 0 ]
    //
    // where theta is chosen to ensure linear non-increase (e1.d_alpha1 +
    // e2.d_alpha2 < 0) and also to make sure that a bound is hit.
    // Specifically:
    //
    // theta = -1.1 * |h1| : if e2 >= e1 and alpha1 is positive
    // theta = -1.1 * |v1| : if e2 >= e1 and alpha1 is negative
    // theta = +1.1 * |h1| : if e2 <  e1 and alpha1 is positive
    // theta = +1.1 * |v1| : if e2 <  e1 and alpha1 is negative
    //

    L_DOUBLE K11;
    L_DOUBLE K22;
    L_DOUBLE K12;
    L_DOUBLE negdetH;
    L_DOUBLE f;

    K11 = problem.opt_get_G_tau(i1-1,i1-1);
    K22 = problem.opt_get_G_tau(i2-1,i2-1);
    K12 = problem.opt_get_G_tau(i1-1,i2-1);

    negdetH = K11+K22-(2.0*K12);

    f = problem.f;

    if ( ( negdetH > ZZERO ) || ( negdetH < -ZZERO ) )
    {
        /*
           Nonsingular case
        */

        #ifdef DEBUG
        std::cerr << "noningular step\n";
        #endif

        if ( f_zero )
        {
            d_b = (((K12-K22)*e1)+((K12-K11)*e2))/negdetH;

            d_alpha1 = (-e1+e2)/negdetH;
            d_alpha2 = (+e1-e2)/negdetH;
        }

        else
        {
            d_b = ((((K11*K22)-(K12*K12))*f)+((K12-K22)*e1)+((K12-K11)*e2))/negdetH;

            d_alpha1 = (((K12-K22)*f)-e1+e2)/negdetH;
            d_alpha2 = (((K12-K11)*f)+e1-e2)/negdetH;
        }

        d_f = -(problem.f);

        #ifdef DEBUG
        std::cerr << "Check me:\n";
        std::cerr << "[ da1 ]     [ K11 K12 1 ]-1 [ e1 ]\n";
        std::cerr << "[ da2 ] = - [ K12 K22 1 ]   [ e2 ]\n";
        std::cerr << "[ db  ]     [ 1   1   0 ]   [ f  ]\n";
        std::cerr << "e1   = " << e1   << "\n";
        std::cerr << "e2   = " << e2   << "\n";
        std::cerr << "f    = " << f    << "\n";
        std::cerr << "K11  = " << K11  << "\n";
        std::cerr << "K12  = " << K12  << "\n";
        std::cerr << "K22  = " << K22  << "\n";
        std::cerr << "da1  = " << d_alpha1 << "\n";
        std::cerr << "da2  = " << d_alpha2 << "\n";
        std::cerr << "db   = " << d_b  << "\n";
        std::cerr << "d_f  = " << d_f  << "\n";
        #endif

        d_J_epart = -((e1-e2)*(e1-e2))/negdetH;
    }

    else
    {
        /*
           Singular case
        */

        #ifdef DEBUG
        std::cerr << "singular step " << e1 << "," << e2 << "\n";
        #endif

        L_DOUBLE theta;

        d_b = K11-K12;

        d_alpha1 = -1.0;
        d_alpha2 = +1.0;

        d_f = 0.0;

        if ( tau2 == +1 )
        {
            if ( e2 >= e1 )
            {
                theta = -1.1 * (problem.fast_h)[i2];
            }

            else
            {
                theta = 1.1 * (problem.fast_h)[i2];
            }
        }

        else
        {
            if ( e2 >= e1 )
            {
                theta = -1.1 * -(problem.fast_v)[i2];
            }

            else
            {
                theta = 1.1 * -(problem.fast_v)[i2];
            }
        }

        d_b *= theta;

        d_alpha1 *= theta;
        d_alpha2 *= theta;

        d_J_epart = (d_alpha1*e1) + (d_alpha2*e2);
    }

    /*
       Now let's scale the step
    */

    L_DOUBLE step_scale = 1.0;
    L_DOUBLE ceta;

    f_zero_next = 1;

    if ( ( d_alpha1 > -ZZERO/10 ) && ( d_alpha1 < ZZERO/10 ) )
    {
        d_alpha1 = 0.0;
    }

    if ( ( d_alpha2 > -ZZERO/10 ) && ( d_alpha2 < ZZERO/10 ) )
    {
        d_alpha2 = 0.0;
    }

    if ( d_alpha1 < 0.0 )
    {
        if ( tau1 > 0 )
        {
            ceta = -((problem.fast_alpha)[i1-1]) / d_alpha1;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }

        else
        {
            ceta = ( (problem.fast_v)[i1-1] - ((problem.fast_alpha)[i1-1]) ) / d_alpha1;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }
    }

    else if ( d_alpha1 > 0.0 )
    {
        if ( tau1 > 0 )
        {
            ceta = ( (problem.fast_h)[i1-1] - ((problem.fast_alpha)[i1-1]) ) / d_alpha1;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }

        else
        {
            ceta = -((problem.fast_alpha)[i1-1]) / d_alpha1;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }
    }

    if ( d_alpha2 < 0.0 )
    {
        if ( tau2 > 0 )
        {
            ceta = -((problem.fast_alpha)[i2-1]) / d_alpha2;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }

        else
        {
            ceta = ( (problem.fast_v)[i2-1] - ((problem.fast_alpha)[i2-1]) ) / d_alpha2;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }
    }

    else if ( d_alpha2 > 0.0 )
    {
        if ( tau2 > 0 )
        {
            ceta = ( (problem.fast_h)[i2-1] - ((problem.fast_alpha)[i2-1]) ) / d_alpha2;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }

        else
        {
            ceta = -((problem.fast_alpha)[i2-1]) / d_alpha2;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }
    }

    #ifdef DEBUG
    std::cerr << "scaling = " << step_scale << "\n";
    #endif

    d_b *= step_scale;

    d_alpha1 *= step_scale;
    d_alpha2 *= step_scale;

    d_f *= step_scale;

    if ( ( negdetH > ZZERO ) || ( negdetH < -ZZERO ) )
    {
        /*
           Step in quadratically curved direction, so scaling...
        */

        d_J_epart *= (step_scale-((step_scale*step_scale)/2.0));
    }

    else
    {
        /*
           Step in linear direction, so scaling is linear.
        */

        d_J_epart *= step_scale;
    }

    /*
       Return zero if step too small.  This is from Platt, but can only
       be applied validly if f is zero (as it assumes that d_alpha1 and
       d_alpha2 are the same up to sign).
    */

    if ( f_zero && ( fabs(d_alpha1) < sol_tol*(fabs((2.0*((problem.fast_alpha)[i1-1]))+d_alpha1)+sol_tol) ) )
    {
        return 0;
    }

    return 1;
}

int actually_take_step_SMO(L_DOUBLE &d_J_epart, long i1, long i2, int tau1, int tau2, L_DOUBLE &d_alpha1, L_DOUBLE &d_alpha2, L_DOUBLE &d_b, L_DOUBLE &d_f, int &f_zero, int &f_zero_next, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J)
{
    #ifdef DEBUG
    std::cerr << "d_b = " << d_b << "\n";
    std::cerr << "d_f = " << d_f << "\n";
    std::cerr << "i1 = " << i1 << "\n";
    std::cerr << "i2 = " << i2 << "\n";
    std::cerr << "d_alpha1 = " << d_alpha1 << "\n";
    std::cerr << "d_alpha2 = " << d_alpha2 << "\n";
    #endif

    /*
       ...but first, we may need to unconstrain alpha1 and/or alpha2
    */

    if ( ( (problem.fast_tau)[i1-1] != +1 ) && ( (problem.fast_tau)[i1-1] != -1 ) )
    {
        if ( tau1 == +1 )
        {
            (problem.fast_tau)[i1-1] = +1;
        }

        else
        {
            (problem.fast_tau)[i1-1] = -1;
        }
    }

    if ( ( (problem.fast_tau)[i2-1] != +1 ) && ( (problem.fast_tau)[i2-1] != -1 ) )
    {
        if ( tau2 == +1 )
        {
            (problem.fast_tau)[i2-1] = +1;
        }

        else
        {
            (problem.fast_tau)[i2-1] = -1;
        }
    }

    #ifdef SMO_DEBUG
    std::cerr << "===============================================\n";
    std::cerr << "before:\n";
    E_CHECKER_QUIET
    #endif

    (*J) += d_J_epart;

    /*problem.opt_step_two(d_b,d_f,d_alpha1,d_alpha2,d_e1,d_e2,i1-1,i2-1);*/
    problem.opt_step_two(d_b,d_f,d_alpha1,d_alpha2,i1-1,i2-1);

    #ifdef SMO_DEBUG
    std::cerr << "after:\n";
    E_CHECKER_QUIET
    #endif   

/*    #ifdef DEBUG
    std::cerr << problem << "\n";
    #endif
*/

    f_zero = f_zero_next;

    /*
       Be careful with constraints.  If alpha wanders too close to the
       boundary, reel it in and constrain it.
    */

    if ( ( (problem.fast_alpha)[i1-1] <= ZZERO ) && ( (problem.fast_alpha)[i1-1] >= -ZZERO ) )
    {
        (problem.fast_alpha)[i1-1] = 0.0;

        problem.opt_constrain_Z(i1-1);
    }

    else if ( ( (problem.fast_alpha)[i1-1] <= ((problem.fast_v)[i1-1])+ZZERO ) && ( (problem.fast_tau)[i1-1] == -1 ) )
    {
        (problem.fast_alpha)[i1-1] = (problem.fast_v)[i1-1];
        (problem.fast_tau)[i1-1]   = -2;
    }

    else if ( ( (problem.fast_alpha)[i1-1] >= ((problem.fast_h)[i1-1])-ZZERO ) && ( (problem.fast_tau)[i1-1] == +1 ) )
    {
        (problem.fast_alpha)[i1-1] = (problem.fast_h)[i1-1];
        (problem.fast_tau)[i1-1]   = +2;
    }

    if ( ( (problem.fast_alpha)[i2-1] <= ZZERO ) && ( (problem.fast_alpha)[i2-1] >= -ZZERO ) )
    {
        (problem.fast_alpha)[i2-1] = 0.0;

        problem.opt_constrain_Z(i2-1);
    }

    else if ( ( (problem.fast_alpha)[i2-1] <= ((problem.fast_v)[i2-1])+ZZERO ) && ( (problem.fast_tau)[i2-1] == -1 ) )
    {
        (problem.fast_alpha)[i2-1] = (problem.fast_v)[i2-1];
        (problem.fast_tau)[i2-1]   = -2;
    }

    else if ( ( (problem.fast_alpha)[i2-1] >= ((problem.fast_h)[i2-1])-ZZERO ) && ( (problem.fast_tau)[i2-1] == +1 ) )
    {
        (problem.fast_alpha)[i2-1] = (problem.fast_h)[i2-1];
        (problem.fast_tau)[i2-1]   = +2;
    }

    return 1;

    /* burn warnings burn */

    sol_tol = 0.0;
}



unsigned long solve_SMO_fixed_bias(SVdata &problem, L_DOUBLE &sol_tol, unsigned long epochs, volatile int *async_exit_flag, volatile unsigned long *loop_count)
{
    long i;
    long numChanged = 0;
    int examineAll = 1;

    //
    // Main algorithmic loop.
    //

    // 
    // The fixed bias version of SMO is essentially trivial.  The active
    // set size is 1, so there is no need for a second choice heuristic
    // and the related complications.  Furthermore, even the singular case
    // becomes trivial to evaluate.
    //

    for ( ; ( ( (*loop_count) <= epochs ) && ( !(*async_exit_flag) ) ) ; )
    {
        numChanged = 0;

        for ( i = 1 ; i <= problem.N ; i++ )
        {
            switch ( (problem.fast_tau)[i-1] )
            {
                case 0:
                {
                    /*
                       point constrained at zero.  If examineAll is set then
                       proceed to examine all posible alternative "states"
                       for this variable.
                    */

                    if ( examineAll )
                    {
                        switch ( (problem.fast_contype)[i-1] )
                        {
                            case 0:
                            {
                                break;
                            }

                            case 1:
                            {
                                numChanged += takeStep_SMO_fixed_bias(i,-1,problem,sol_tol);

                                break;
                            }

                            case 2:
                            {
                                numChanged += takeStep_SMO_fixed_bias(i,+1,problem,sol_tol);

                                break;
                            }

                            case 3:
                            {
                                if ( rand() % 2 )
                                {
                                    if ( takeStep_SMO_fixed_bias(i,-1,problem,sol_tol) )
                                    {
                                        numChanged++;
                                    }

                                    else
                                    {
                                        numChanged += takeStep_SMO_fixed_bias(i,+1,problem,sol_tol);
                                    }
                                }

                                else
                                {
                                    if ( takeStep_SMO_fixed_bias(i,+1,problem,sol_tol) )
                                    {
                                        numChanged++;
                                    }

                                    else
                                    {
                                        numChanged += takeStep_SMO_fixed_bias(i,-1,problem,sol_tol);
                                    }
                                }

                                break;
                            }

                            default:
                            {
                                L_THROW(0);

                                break;
                            }
                        }
                    }

                    break;
                }

                case -1:
                {
                    /*
                       Variable is free, so examine step
                    */

                    numChanged += takeStep_SMO_fixed_bias(i,-1,problem,sol_tol);

                    break;
                }

                case +1:
                {
                    /*
                       Variable is free, so examine step
                    */

                    numChanged += takeStep_SMO_fixed_bias(i,+1,problem,sol_tol);

                    break;
                }

                case +2:
                {
                    /*
                       point constrained at h.  If examineAll is set then
                       proceed to examine freeing this variable.
                    */

                    if ( examineAll )
                    {
                        numChanged += takeStep_SMO_fixed_bias(i,+1,problem,sol_tol);
                    }

                    break;
                }

                case -2:
                {
                    /*
                       point constrained at v.  If examineAll is set then
                       proceed to examine freeing this variable.
                    */

                    if ( examineAll )
                    {
                        numChanged += takeStep_SMO_fixed_bias(i,-1,problem,sol_tol);
                    }

                    break;
                }

                default:
                {
                    L_THROW(0);

                    break;
                }
            }
        }

        (*loop_count) += numChanged;

        if ( examineAll )
        {
            if ( numChanged == 0 )
            {
                goto getout;
            }

            examineAll = 0;
        }

        else if ( numChanged == 0 )
        {
            examineAll = 1;
        }
    }

    getout:

    return (*loop_count);
}

int takeStep_SMO_fixed_bias(long i, int tau, SVdata &problem, L_DOUBLE &sol_tol)
{
    L_DOUBLE alpha;
    L_DOUBLE e;
    L_DOUBLE vbar;
    L_DOUBLE hbar;

    alpha = (problem.fast_alpha)[i-1];
    e     = ( tau == +1 ) ? problem.opt_get_e_pos(i-1)         : problem.opt_get_e_neg(i-1);
    vbar  = ( tau == +1 ) ? 0.0                                : ( (double) (problem.fast_v)[i-1] );
    hbar  = ( tau == +1 ) ? ( (double) (problem.fast_h)[i-1] ) : 0.0;

    if ( ( ( e < -sol_tol ) && ( alpha < hbar ) ) || ( ( e > sol_tol ) && ( alpha > vbar ) ) )
    {
        L_DOUBLE K11;
        L_DOUBLE d_alpha;
        L_DOUBLE d_e;

        K11 = problem.opt_get_G_tau(i-1,i-1);

        if ( ( K11 > ZZERO ) || ( K11 < -ZZERO ) )
        {
            /*
               Nonsingular case.  This is just a newton step.
            */

            d_alpha = -e/K11;
            d_e     = -e;
        }

        else
        {
            /*
               Singular case.  Choose linear descent and head to boundary.
            */

            if ( e < 0.0 )
            {
                d_alpha = (hbar-vbar);
            }

            else
            {
                d_alpha = -(hbar-vbar);
            }

            d_e = 0.0;
        }

        /*
           Now let's scale the step
        */

        L_DOUBLE step_scale = 1.0;
        L_DOUBLE ceta;

        if ( ( d_alpha > -ZZERO/10 ) && ( d_alpha < ZZERO/10 ) )
        {
            d_alpha = 0.0;
        }

        if ( d_alpha < 0.0 )
        {
            ceta = ( vbar - alpha ) / d_alpha;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
            }
        }

        else if ( d_alpha > 0.0 )
        {
            ceta = ( hbar - alpha ) / d_alpha;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
            }
        }

        d_alpha *= step_scale;
        d_e *= step_scale;

        /*
           Return zero if step too small.  This is from Platt, approximately.
        */

        if ( fabs(d_alpha) < sol_tol*(fabs((2.0*alpha)+d_alpha)+sol_tol) )
        {
            return 0;
        }

        /*
           Unconstrain if necessary
        */

        if ( ( (problem.fast_tau)[i-1] != +1 ) && ( (problem.fast_tau)[i-1] != -1 ) )
        {
            if ( tau == +1 )
            {
                (problem.fast_tau)[i-1] = +1;
            }

            else
            {
                (problem.fast_tau)[i-1] = -1;
            }
        }

        /*
           Take the step;
        */

/*        problem.opt_step_one(d_alpha,d_e,i-1);*/
        problem.opt_step_one(d_alpha,i-1);

        /*
           Be careful with constraints.  If alpha wanders too close to the
           boundary, real it in and constrain it.
        */

        if ( ( (problem.fast_alpha)[i-1] <= ZZERO ) && ( (problem.fast_alpha)[i-1] >= -ZZERO ) )
        {
            (problem.fast_alpha)[i-1] = 0.0;

            problem.opt_constrain_Z(i-1);
        }

        else if ( ( (problem.fast_alpha)[i-1] <= ((problem.fast_v)[i-1])+ZZERO ) && ( (problem.fast_tau)[i-1] == -1 ) )
        {
            (problem.fast_alpha)[i-1] = (problem.fast_v)[i-1];
            (problem.fast_tau)[i-1]   = -2;
        }

        else if ( ( (problem.fast_alpha)[i-1] >= ((problem.fast_h)[i-1])-ZZERO ) && ( (problem.fast_tau)[i-1] == +1 ) )
        {
            (problem.fast_alpha)[i-1] = (problem.fast_h)[i-1];
            (problem.fast_tau)[i-1]   = +2;
        }

        return 1;
    }

    return 0;
}






/***********************************************************************

                            D2C optimiser
                            =============

***********************************************************************/

/*
   Size of buffer used when sorting the most non-optimal points.
*/

/*
   SORTINTOPOS(ee,new_tau,alpha_dist):

   - ee is the gradient at alpha
   - new_tau is what tau will become if we move in the -ve alpha direction.
   - alpha_dist is the distance we can move in the -ve alpha direction
     before hitting a bound.

   This macro will attempt to put variable i into the buffer of maximal
   gradients in the -ve alpha direction - that is, a buffer ordered in
   terms of decreasing gradients, starting with the most positive gradient.


   SORTINTONEG(ee,new_tau,alpha_dist):

   - ee is the gradient at alpha
   - new_tau is what tau will become if we move in the +ve alpha direction.
   - alpha_dist is the distance we can move in the +ve alpha direction
     before hitting a bound.

   This macro will attempt to put variable i into the buffer of maximal
   gradients in the +ve alpha direction - that is, a buffer ordered in
   terms of increasing gradients, starting with the most negative gradient.
*/

#define MACROREM(a)

#define SORTINTOPOS(_ee_,_old_tau_,_new_tau_)                           \
{                                                                       \
    j = problem.d2c_buff_size;                                          \
                                                                        \
    MACROREM("The second condition here overrides the first cond");     \
    MACROREM("if the position in the buffer is empty.");                \
                                                                        \
    while ( ( ( _ee_ > max_pos_val[j-1] ) || !(max_pos_pos[j-1]) ) && j ) \
    {                                                                   \
        j--;                                                            \
    }                                                                   \
                                                                        \
    if ( j < problem.d2c_buff_size )                                    \
    {                                                                   \
        k = (problem.d2c_buff_size)-1;                                  \
                                                                        \
        while ( k > j )                                                 \
        {                                                               \
            max_pos_pos[k]        = max_pos_pos[k-1];                   \
            max_pos_val[k]        = max_pos_val[k-1];                   \
            max_pos_old_tau[k]    = max_pos_old_tau[k-1];               \
            max_pos_new_tau[k]    = max_pos_new_tau[k-1];               \
                                                                        \
            k--;                                                        \
        }                                                               \
                                                                        \
        max_pos_pos[j]        = i;                                      \
        max_pos_val[j]        = _ee_;                                   \
        max_pos_old_tau[j]    = _old_tau_;                              \
        max_pos_new_tau[j]    = _new_tau_;                              \
                                                                        \
        if ( pos_buffer_len < problem.d2c_buff_size )                   \
        {                                                               \
            pos_buffer_len++;                                           \
        }                                                               \
    }                                                                   \
}

#define SORTINTONEG(_ee_,_old_tau_,_new_tau_)                           \
{                                                                       \
    j = problem.d2c_buff_size;                                          \
                                                                        \
    while ( ( ( _ee_ < max_neg_val[j-1] ) || !(max_neg_pos[j-1]) ) && j ) \
    {                                                                   \
        j--;                                                            \
    }                                                                   \
                                                                        \
    if ( j < problem.d2c_buff_size )                                    \
    {                                                                   \
        k = (problem.d2c_buff_size)-1;                                  \
                                                                        \
        while ( k > j )                                                 \
        {                                                               \
            max_neg_pos[k]        = max_neg_pos[k-1];                   \
            max_neg_val[k]        = max_neg_val[k-1];                   \
            max_neg_old_tau[k]    = max_neg_old_tau[k-1];               \
            max_neg_new_tau[k]    = max_neg_new_tau[k-1];               \
                                                                        \
            k--;                                                        \
        }                                                               \
                                                                        \
        max_neg_pos[j]        = i;                                      \
        max_neg_val[j]        = _ee_;                                   \
        max_neg_old_tau[j]    = _old_tau_;                              \
        max_neg_new_tau[j]    = _new_tau_;                              \
                                                                        \
        if ( neg_buffer_len < problem.d2c_buff_size )                   \
        {                                                               \
            neg_buffer_len++;                                           \
        }                                                               \
    }                                                                   \
}

/*
               Heuristic 1 - if largest violator at bound, find largest
                             complimentary violator at bound and step in
                             both.

                             QUESTION: This is assuming a consistent upper
                                       (lower) bound on alpha.  However, this
                                       won't necessarily hold - what's the
                                       best way to deal with this?  What is
                                       the best method of choosing in this
                                       case?

                             ANSWER?: if the distance alpha_i may travel is
                                      a_i then the change in the objective fn
                                      will be (up to a negative scale
                                      factor):
                                      a_i.e_i + a_j.e_j.
                                      We want to maximise this magnitude, and
                                      to do so we maximise a_i.e_i and
                                      a_j.e_j.  Of course, need to ensure
                                      that the directions of the steps are
                                      complementary - ie.
                                      - if alpha_i = v => either alpha_j = h
                                        or alpha_j = 0 and tau_j = -1
                                      - if alpha_i = h => either alpha_j = v
                                        or alpha_j = 0 and tau_j = +1
                                      - if alpha_i = 0 and tau_i = +1 then
                                        either alpha_j = h or alpha_j = 0
                                        and tau_j = -1
                                      - if alpha_i = 0 and tau_i = -1 then
                                        either alpha_j = v or alpha_j = 0
                                        and tau_j = +1

                             OBVIOUS FOLLOWUP Q: given that we want to max
                                                 a_i.e_i and a_j.e_j, why not
                                                 just do this straight off?
                                                 To put it another way - am I
                                                 looking at the correct
                                                 metric here when ordering
                                                 the e_i values?  Would it
                                                 be better to order the
                                                 a_i.e_i values instead?

                             FOR NOW: just blindly follow the algorithm,
                                      ignore the magnitude of a_i.

                             IMPORTANT: when f != 0 then it is necessary to
                                        check to ensure that a step has
                                        actually occured, or infinite loops
                                        may result.


               Heuristic 2 - if largest violator not at bound, do an
                             exhaustive search to find the largest change in
                             the objective (ASSUMPTION - ignore f component,
                             base only on e) and take this step.

               Heuristic 3 - the first 2 heuristics can fail if f != 0.  In
                             either case, simply revert to SMO approximately
                             when examineall = 1.  This is NOT the optimal
                             thing to do.
*/


#define E_CHECKER_GEN                                                   \
{                                                                       \
    L_DOUBLE __new_e;                                                   \
    long __I;                                                           \
    long __J;                                                           \
    std::cerr << "phantom - echeck\n";                                  \
    for ( __I = 1 ; __I <= problem.N ; __I++ )                          \
    {                                                                   \
        __new_e = problem.b - (problem.fast_z)[__I-1];                  \
                                                                        \
        for ( __J = 1 ; __J <= problem.N ; __J++ )                      \
        {                                                               \
            if ( problem.fast_tau[__J-1] != 0 )                         \
            {                                                           \
                __new_e += (problem.opt_get_G_tau(__I-1,__J-1)) * (problem.fast_alpha)[__J-1]; \
            }                                                           \
        }                                                               \
                                                                        \
        if ( (problem.fast_contype)[__I-1] == 2 )                       \
        {                                                               \
            __new_e += (problem.fast_rho)[__I-1];                       \
        }                                                               \
                                                                        \
        if ( (problem.fast_contype)[__I-1] == 1 )                       \
        {                                                               \
            __new_e -= (problem.fast_rho_star)[__I-1];                  \
        }                                                               \
                                                                        \
        if ( fabs((problem.fast_e_tau[__I-1])-__new_e) > 1e-6 )         \
        {                                                               \
            std::cerr << "fuck fuck fuck fuck\n";                       \
            exit(1);                                                    \
        }                                                               \
    }                                                                   \
}


unsigned long solve_d2c(SVdata &problem, L_DOUBLE &sol_tol, unsigned long epochs, volatile int *async_exit_flag, volatile unsigned long *loop_count, JTYPE *J)
{
  int step_taken;
  L_DOUBLE ee;
  L_DOUBLE bj_da_dJ = 0.0;
  L_DOUBLE d_J_epart = 0.0;
  int f_zero;
  int f_zero_next = 0;
  int first_try;
  long i,j,k,l,m;
  long i1,i2 = 0;
  int oldtau1,oldtau2;
  int newtau1,newtau2 = 0;
  L_DOUBLE e1,e2;
  long pos_buffer_len;
  long neg_buffer_len;
  long     *max_pos_pos;
  L_DOUBLE *max_pos_val;
  int      *max_pos_old_tau;
  int      *max_pos_new_tau;
  long     *max_neg_pos;
  L_DOUBLE *max_neg_val;
  int      *max_neg_old_tau;
  int      *max_neg_new_tau;
  long     *max_sec_pos;
  L_DOUBLE *max_sec_val;
  int      *max_sec_old_tau;
  int      *max_sec_new_tau;
  L_DOUBLE d_alpha1 = 0.0;
  L_DOUBLE d_alpha2 = 0.0;
  L_DOUBLE d_b = 0.0;
  L_DOUBLE d_f = 0.0;
  L_DOUBLE d_alpha1_trial = 0.0;
  L_DOUBLE d_alpha2_trial = 0.0;
  L_DOUBLE d_b_trial = 0.0;
  L_DOUBLE d_f_trial = 0.0;

  REPNEWB(max_pos_pos,long,(problem.d2c_buff_size));
  REPNEWB(max_pos_val,L_DOUBLE,(problem.d2c_buff_size));
  REPNEWB(max_pos_old_tau,int,(problem.d2c_buff_size));
  REPNEWB(max_pos_new_tau,int,(problem.d2c_buff_size));
  REPNEWB(max_neg_pos,long,(problem.d2c_buff_size));
  REPNEWB(max_neg_val,L_DOUBLE,(problem.d2c_buff_size));
  REPNEWB(max_neg_old_tau,int,(problem.d2c_buff_size));
  REPNEWB(max_neg_new_tau,int,(problem.d2c_buff_size));

  f_zero = 0;

  if ( ( problem.f > -sol_tol ) && ( problem.f < sol_tol ) )
  {
    f_zero = 1;
  }

  if ( !D2CSMO_OPTIM(problem.svflags) )
  {
    /* SEE ALSO OPTIMISED VERSION BELOW */
    /*
       max_pos_pos: descent in the +alpha direction
       max_neg_pos: descent in the -alpha direction

       *_new_tau: tau value used to calculate this value.
    */

    for ( ; ( ( (*loop_count) <= epochs ) && ( !(*async_exit_flag) ) ) ; (*loop_count)++ )
    {
        /*
           First thing - sort the most non-optimal into positive and
           negative buffers.
        */

        #ifdef DEBUG
        std::cerr << "Clearing buffers\n";
        #endif

        for ( i = 1 ; i <= problem.d2c_buff_size ; i++ )
        {
            max_pos_pos[i-1] = 0; /* 0 indicates unfilled */
            max_neg_pos[i-1] = 0; /* 0 indicates unfilled */

            max_pos_val[i-1] = 0.0;
            max_neg_val[i-1] = 0.0;
        }

        #ifdef DEBUG
        std::cerr << "Filling buffers\n";
        #endif

        pos_buffer_len = 0;
        neg_buffer_len = 0;

        for ( i = 1 ; i <= problem.N ; i++ )
        {
            switch ( (problem.fast_tau)[i-1] )
            {
                /*
                    tau | contype | alpha can? | gradients
                   -----+---------+------------+-----------
                     +2 | na      | decrease   | +ve
                     -2 | na      | increase   | -ve
                     +1 | na      | either     | +ve/-ve
                     -1 | na      | either     | +ve/-ve
                     0  | 0       | nothing    | na
                        | 1       | decrease   | +ve
                        | 2       | increase   | -ve
                        | 3       | either     | +ve/-ve
                */

                case +2: { ee = problem.opt_get_e_pos(i-1); SORTINTOPOS(ee,+2,+1); break; }
                case -2: { ee = problem.opt_get_e_neg(i-1); SORTINTONEG(ee,-2,-1); break; }
                case +1: { ee = problem.opt_get_e_pos(i-1); SORTINTOPOS(ee,+1,+1);
                                                            SORTINTONEG(ee,+1,+1); break; }
                case -1: { ee = problem.opt_get_e_neg(i-1); SORTINTOPOS(ee,-1,-1);
                                                            SORTINTONEG(ee,-1,-1); break; }

                case 0:
                {
                    switch ( (problem.fast_contype)[i-1] )
                    {
                        case 0: { break; }
                        case 1: { ee = problem.opt_get_e_neg(i-1); SORTINTOPOS(ee,0,-1); break; }
                        case 2: { ee = problem.opt_get_e_pos(i-1); SORTINTONEG(ee,0,+1); break; }
                        case 3: { ee = problem.opt_get_e_neg(i-1); SORTINTOPOS(ee,0,-1);
                                  ee = problem.opt_get_e_pos(i-1); SORTINTONEG(ee,0,+1); break; }

                        default:
                        {
                            L_THROW(0);

                            break;
                        }
                    }

                    break;
                }

                default:
                {
                    L_THROW(0);

                    break;
                }
            }
        }

        /*
           Is the solution optimal?  This checks the violation gap and also
           (implicitly) the optimality of points not at the boundary.  Also,
           need to check that f is zero here.
        */

        #ifdef DEBUG
        std::cerr << "Testing optimality\n";
        #endif

        if ( f_zero && ( !(max_pos_pos[0]) || ( max_pos_val[0] <  sol_tol ) )
                    && ( !(max_neg_pos[0]) || ( max_neg_val[0] > -sol_tol ) ) )
        {
            goto exit_point;
        }

        step_taken = 0;

        i = 0;
        j = 0;

        /*
           Heuristics 1 and 2 will be used until either a step has been
           taken or the lists themselves are both exhausted.
        */

        #ifdef DEBUG
        std::cerr << "Beginning heuristic 1/2 loop...\n";
        #endif

        while ( ( !step_taken        ) &&
                ( i < pos_buffer_len ) &&
                ( j < neg_buffer_len )    )
        {
            #ifdef DEBUG
            std::cerr << "Finding startpoint\n";
            #endif

            if ( max_pos_val[i] > -max_neg_val[j] )
            {
                i1      = max_pos_pos[i];
                e1      = max_pos_val[i];
                oldtau1 = max_pos_old_tau[i];
                newtau1 = max_pos_new_tau[i];

                i++;

                max_sec_pos        = max_neg_pos;
                max_sec_val        = max_neg_val;
                max_sec_new_tau    = max_neg_new_tau;
                max_sec_old_tau    = max_neg_new_tau;

                k = j;
                l = neg_buffer_len;
            }

            else
            {
                i1      = max_neg_pos[j];
                e1      = max_neg_val[j];
                oldtau1 = max_neg_old_tau[j];
                newtau1 = max_neg_new_tau[j];

                j++;

                max_sec_pos        = max_pos_pos;
                max_sec_val        = max_pos_val;
                max_sec_new_tau    = max_pos_new_tau;
                max_sec_old_tau    = max_pos_new_tau;

                k = i;
                l = pos_buffer_len;
            }

            /*
               Heuristic 1
            */

            if ( ( oldtau1 == -2 ) ||
                 ( oldtau1 ==  0 ) ||
                 ( oldtau1 == +2 )    )
            {
                #ifdef DEBUG
                std::cerr << "Heuristic 1\n";
                #endif

                /*
                   Scan through the secondary list.  The first to be found
                   with oldtau begin one of -2,0,+2 is at a bound and
                   pointed in the correct direction.  Hence attempt to take
                   a step w.r.t. this one.  If this attempted step succeeds
                   then everything is good, we can exit the while loop.  If
                   this step fails, keep looking until another potential is
                   found or fall-through occurs.
                */

                for ( m = k ; ( ( m < l ) && !step_taken ) ; m++ )
                {
                    if ( ( max_sec_old_tau[m] == -2 ) ||
                         ( max_sec_old_tau[m] ==  0 ) ||
                         ( max_sec_old_tau[m] == +2 )    )
                    {
                        step_taken = takeStep_d2c(i1,max_sec_pos[m],newtau1,max_sec_new_tau[m],e1,max_sec_val[m],f_zero,problem,sol_tol,J);

                        #ifdef DEBUG
                        if ( step_taken )
                        {
                            std::cerr << "**** Actually took step " << i1 << " - " << max_sec_pos[m] << "\n";
                        }
                        #endif
                    }
                }
            }

            /*
               Heuristic 2
            */

            if ( !step_taken )
            {
                #ifdef DEBUG
                std::cerr << "Heuristic 2\n";
                #endif

                /*
                   Find the largest feasible step by searching through the
                   secondary gradient list.
                */

                first_try = 1;

                for ( m = k ; m < l ; m++ )
                {
                    /*
                       Trial the step - this will record the step itself, and
                       also calculate the change in the objective function J
                       due to the step (it will return 0 if no step is
                       possible).
                    */

                    if ( i1 != max_sec_pos[m] )
                    {
                        if ( trial_step_d2c(d_J_epart,i1,max_sec_pos[m],newtau1,max_sec_new_tau[m],e1,max_sec_val[m],d_alpha1_trial,d_alpha2_trial,d_b_trial,d_f_trial,f_zero,f_zero_next,problem,sol_tol) )
                        {
                            /*
                               OK, the step is feasible.  If either this is the
                               first feasible step (first_try = 1) or it results
                               in the largest objective change then record it.
                               Otherwise, ignore it.
                            */

                            if ( first_try )
                            {
                                bj_da_dJ = d_J_epart;

                                i2      = max_sec_pos[m];
                                e2      = max_sec_val[m];
                                oldtau2 = max_sec_old_tau[m];
                                newtau2 = max_sec_new_tau[m];

                                d_alpha1 = d_alpha1_trial;
                                d_alpha2 = d_alpha2_trial;
                                d_b      = d_b_trial;
                                d_f      = d_f_trial;

                                first_try = 0;
                            }

                            else if ( d_J_epart < bj_da_dJ )
                            {
                                bj_da_dJ = d_J_epart;

                                i2      = max_sec_pos[m];
                                e2      = max_sec_val[m];
                                oldtau2 = max_sec_old_tau[m];
                                newtau2 = max_sec_new_tau[m];

                                d_alpha1 = d_alpha1_trial;
                                d_alpha2 = d_alpha2_trial;
                                d_b      = d_b_trial;
                                d_f      = d_f_trial;
                            }
                        }
                    }
                }

                /*
                   Take the largest feasible step, if there is one.  The
                   function actually_take_step_d2c will always return 1,
                   so this will correctly set step_taken.
                */

                if ( !first_try )
                {
                    #ifdef DEBUG
                    std::cerr << "**** Actually take step " << i1 << " - " << i2 << "\n";
                    #endif

                    step_taken = actually_take_step_d2c(bj_da_dJ,i1,i2,newtau1,newtau2,d_alpha1,d_alpha2,d_b,d_f,f_zero,f_zero_next,problem,sol_tol,J);
                }
            }
        }

        if ( !step_taken )
        {
            #ifdef DEBUG
            std::cerr << "Heuristic 3 loop\n";
            #endif

            /*
               Heuristic 3
            */

            for ( m = 1 ; ( ( m <= problem.N ) && !step_taken ) ; m++ )
            {
                switch ( (problem.fast_tau)[m-1] )
                {
                    case 0:
                    {
                        switch ( (problem.fast_contype)[m-1] )
                        {
                            case 0: { break; }

                            case 1: { step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J); break; }
                            case 2: { step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J); break; }

                            case 3:
                            {
                                if ( rand() % 2 )
                                {
                                    if ( !( step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J) ) )
                                    {
                                        step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J);
                                    }
                                }

                                else
                                {
                                    if ( !( step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J) ) )
                                    {
                                        step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J);
                                    }
                                }

                                break;
                            }

                            default:
                            {
                                L_THROW(0);

                                break;
                            }
                        }

                        break;
                    }

                    case -1: { step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J); break; }
                    case +1: { step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J); break; }
                    case +2: { step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J); break; }
                    case -2: { step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J); break; }
 
                    default:
                    {
                        L_THROW(0);

                        break;
                    }
                }
            }
        }

        /*
           Fall-back case (numerical weirdness has happened).
        */

        if ( !step_taken )
        {
            goto exit_point;
        }

/*        THROW_ASSERT(step_taken);*/
    }
  }

  else if ( !SVM_HINT_PATTERN(problem.svflags) )
  {
    for ( ; ( ( (*loop_count) <= epochs ) && ( !(*async_exit_flag) ) ) ; (*loop_count)++ )
    {
        for ( i = 1 ; i <= problem.d2c_buff_size ; i++ )
        {
            max_pos_pos[i-1] = 0;
            max_neg_pos[i-1] = 0;

            max_pos_val[i-1] = 0.0;
            max_neg_val[i-1] = 0.0;
        }

        pos_buffer_len = 0;
        neg_buffer_len = 0;

        for ( i = 1 ; i <= problem.N ; i++ )
        {
            switch ( (problem.fast_tau)[i-1] )
            {
                case +2: { ee = problem.fast_e_tau[i-1]; SORTINTOPOS(ee,+2,+1); break; }
                case -2: { ee = problem.fast_e_tau[i-1]; SORTINTONEG(ee,-2,-1); break; }
                case +1: { ee = problem.fast_e_tau[i-1]; SORTINTOPOS(ee,+1,+1);
                                                         SORTINTONEG(ee,+1,+1); break; }
                case -1: { ee = problem.fast_e_tau[i-1]; SORTINTOPOS(ee,-1,-1);
                                                         SORTINTONEG(ee,-1,-1); break; }

                default:
                {
                    switch ( (problem.fast_contype)[i-1] )
                    {
                        case 0:  { break; }
                        case 1:  { ee = problem.opt_get_e_neg(i-1); SORTINTOPOS(ee,0,-1); break; }
                        case 2:  { ee = problem.opt_get_e_pos(i-1); SORTINTONEG(ee,0,+1); break; }
                        default: { ee = problem.opt_get_e_neg(i-1); SORTINTOPOS(ee,0,-1);
                                   ee = problem.opt_get_e_pos(i-1); SORTINTONEG(ee,0,+1); break; }
                    }

                    break;
                }
            }
        }

        if ( f_zero && ( !(max_pos_pos[0]) || ( max_pos_val[0] <  sol_tol ) )
                    && ( !(max_neg_pos[0]) || ( max_neg_val[0] > -sol_tol ) ) )
        {
            goto exit_point;
        }

        step_taken = 0;

        i = 0;
        j = 0;

        while ( ( !step_taken        ) &&
                ( i < pos_buffer_len ) &&
                ( j < neg_buffer_len )    )
        {
            if ( max_pos_val[i] > -max_neg_val[j] )
            {
                i1      = max_pos_pos[i];
                e1      = max_pos_val[i];
                oldtau1 = max_pos_old_tau[i];
                newtau1 = max_pos_new_tau[i];

                i++;

                max_sec_pos        = max_neg_pos;
                max_sec_val        = max_neg_val;
                max_sec_new_tau    = max_neg_new_tau;
                max_sec_old_tau    = max_neg_new_tau;

                k = j;
                l = neg_buffer_len;
            }

            else
            {
                i1      = max_neg_pos[j];
                e1      = max_neg_val[j];
                oldtau1 = max_neg_old_tau[j];
                newtau1 = max_neg_new_tau[j];

                j++;

                max_sec_pos        = max_pos_pos;
                max_sec_val        = max_pos_val;
                max_sec_new_tau    = max_pos_new_tau;
                max_sec_old_tau    = max_pos_new_tau;

                k = i;
                l = pos_buffer_len;
            }

            if ( ( oldtau1 == -2 ) ||
                 ( oldtau1 ==  0 ) ||
                 ( oldtau1 == +2 )    )
            {
                for ( m = k ; ( ( m < l ) && !step_taken ) ; m++ )
                {
                    if ( ( max_sec_old_tau[m] == -2 ) ||
                         ( max_sec_old_tau[m] ==  0 ) ||
                         ( max_sec_old_tau[m] == +2 )    )
                    {
                        step_taken = takeStep_d2c(i1,max_sec_pos[m],newtau1,max_sec_new_tau[m],e1,max_sec_val[m],f_zero,problem,sol_tol,J);
                    }
                }
            }

            if ( !step_taken )
            {
                first_try = 1;

                for ( m = k ; m < l ; m++ )
                {
                    if ( i1 != max_sec_pos[m] )
                    {
                        if ( trial_step_d2c(d_J_epart,i1,max_sec_pos[m],newtau1,max_sec_new_tau[m],e1,max_sec_val[m],d_alpha1_trial,d_alpha2_trial,d_b_trial,d_f_trial,f_zero,f_zero_next,problem,sol_tol) )
                        {
                            if ( first_try )
                            {
                                bj_da_dJ = d_J_epart;

                                i2      = max_sec_pos[m];
                                e2      = max_sec_val[m];
                                oldtau2 = max_sec_old_tau[m];
                                newtau2 = max_sec_new_tau[m];

                                d_alpha1 = d_alpha1_trial;
                                d_alpha2 = d_alpha2_trial;
                                d_b      = d_b_trial;
                                d_f      = d_f_trial;

                                first_try = 0;
                            }

                            else if ( d_J_epart < bj_da_dJ )
                            {
                                bj_da_dJ = d_J_epart;

                                i2      = max_sec_pos[m];
                                e2      = max_sec_val[m];
                                oldtau2 = max_sec_old_tau[m];
                                newtau2 = max_sec_new_tau[m];

                                d_alpha1 = d_alpha1_trial;
                                d_alpha2 = d_alpha2_trial;
                                d_b      = d_b_trial;
                                d_f      = d_f_trial;
                            }
                        }
                    }
                }

                if ( !first_try )
                {
                    step_taken = actually_take_step_d2c(bj_da_dJ,i1,i2,newtau1,newtau2,d_alpha1,d_alpha2,d_b,d_f,f_zero,f_zero_next,problem,sol_tol,J);
                }
            }
        }

        if ( !step_taken )
        {
            for ( m = 1 ; ( ( m <= problem.N ) && !step_taken ) ; m++ )
            {
                switch ( (problem.fast_tau)[m-1] )
                {
                    case 0:
                    {
                        switch ( (problem.fast_contype)[m-1] )
                        {
                            case 0: { break; }
                            case 1: { step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J); break; }
                            case 2: { step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J); break; }

                            default:
                            {
                                if ( rand() % 2 )
                                {
                                    if ( !( step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J) ) )
                                    {
                                        step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J);
                                    }
                                }

                                else
                                {
                                    if ( !( step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J) ) )
                                    {
                                        step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J);
                                    }
                                }

                                break;
                            }
                        }

                        break;
                    }

                    case -1: { step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J); break; }
                    case +1: { step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J); break; }
                    case +2: { step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J); break; }
                    default: { step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J); break; }
                }
            }
        }

        if ( !step_taken )
        {
            goto exit_point;
        }
    }
  }

  else if ( (problem.fast_kernel_fn) == NULL )
  {
    #ifdef DEBUG
    E_CHECKER_GEN
    #endif

    for ( ; ( ( (*loop_count) <= epochs ) && ( !(*async_exit_flag) ) ) ; (*loop_count)++ )
    {
        for ( i = 1 ; i <= problem.d2c_buff_size ; i++ )
        {
            max_pos_pos[i-1] = 0;
            max_neg_pos[i-1] = 0;

            max_pos_val[i-1] = 0.0;
            max_neg_val[i-1] = 0.0;
        }

        pos_buffer_len = 0;
        neg_buffer_len = 0;

        for ( i = 1 ; i <= problem.N ; i++ )
        {
            switch ( (problem.fast_tau)[i-1] )
            {
                case +2: { ee = problem.fast_e_tau[i-1]; SORTINTOPOS(ee,+2,+1); break; }
                case -2: { ee = problem.fast_e_tau[i-1]; SORTINTONEG(ee,-2,-1); break; }
                case +1: { ee = problem.fast_e_tau[i-1]; SORTINTOPOS(ee,+1,+1);
                                                         SORTINTONEG(ee,+1,+1); break; }
                case -1: { ee = problem.fast_e_tau[i-1]; SORTINTOPOS(ee,-1,-1);
                                                         SORTINTONEG(ee,-1,-1); break; }

                default:
                {
                    switch ( (problem.fast_contype)[i-1] )
                    {
                        case 0:  { break; }
                        case 1:  { ee = problem.fast_e_tau[i-1]; SORTINTOPOS(ee,0,-1); break; }
                        case 2:  { ee = problem.fast_e_tau[i-1]; SORTINTONEG(ee,0,+1); break; }
                        default: { ee = problem.fast_e_tau[i-1]; SORTINTOPOS(ee,0,-1);
                                                                 SORTINTONEG(ee,0,+1); break; }
                    }

                    break;
                }
            }
        }

        if ( f_zero && ( !(max_pos_pos[0]) || ( max_pos_val[0] <  sol_tol ) )
                    && ( !(max_neg_pos[0]) || ( max_neg_val[0] > -sol_tol ) ) )
        {
            goto exit_point;
        }

        step_taken = 0;

        i = 0;
        j = 0;

        while ( ( !step_taken        ) &&
                ( i < pos_buffer_len ) &&
                ( j < neg_buffer_len )    )
        {
            if ( max_pos_val[i] > -max_neg_val[j] )
            {
                i1      = max_pos_pos[i];
                e1      = max_pos_val[i];
                oldtau1 = max_pos_old_tau[i];
                newtau1 = max_pos_new_tau[i];

                i++;

                max_sec_pos        = max_neg_pos;
                max_sec_val        = max_neg_val;
                max_sec_new_tau    = max_neg_new_tau;
                max_sec_old_tau    = max_neg_new_tau;

                k = j;
                l = neg_buffer_len;
            }

            else
            {
                i1      = max_neg_pos[j];
                e1      = max_neg_val[j];
                oldtau1 = max_neg_old_tau[j];
                newtau1 = max_neg_new_tau[j];

                j++;

                max_sec_pos        = max_pos_pos;
                max_sec_val        = max_pos_val;
                max_sec_new_tau    = max_pos_new_tau;
                max_sec_old_tau    = max_pos_new_tau;

                k = i;
                l = pos_buffer_len;
            }

            if ( ( oldtau1 == -2 ) ||
                 ( oldtau1 ==  0 ) ||
                 ( oldtau1 == +2 )    )
            {
                for ( m = k ; ( ( m < l ) && !step_taken ) ; m++ )
                {
                    if ( ( max_sec_old_tau[m] == -2 ) ||
                         ( max_sec_old_tau[m] ==  0 ) ||
                         ( max_sec_old_tau[m] == +2 )    )
                    {
                        step_taken = takeStep_d2c_pattern(i1,max_sec_pos[m],newtau1,max_sec_new_tau[m],e1,max_sec_val[m],f_zero,problem,sol_tol,J);
                    }
                }
            }

            if ( !step_taken )
            {
                first_try = 1;

                for ( m = k ; m < l ; m++ )
                {
                    if ( i1 != max_sec_pos[m] )
                    {
                        if ( trial_step_d2c(d_J_epart,i1,max_sec_pos[m],newtau1,max_sec_new_tau[m],e1,max_sec_val[m],d_alpha1_trial,d_alpha2_trial,d_b_trial,d_f_trial,f_zero,f_zero_next,problem,sol_tol) )
                        {
                            if ( first_try )
                            {
                                bj_da_dJ = d_J_epart;

                                i2      = max_sec_pos[m];
                                e2      = max_sec_val[m];
                                oldtau2 = max_sec_old_tau[m];
                                newtau2 = max_sec_new_tau[m];

                                d_alpha1 = d_alpha1_trial;
                                d_alpha2 = d_alpha2_trial;
                                d_b      = d_b_trial;
                                d_f      = d_f_trial;

                                first_try = 0;
                            }

                            else if ( d_J_epart < bj_da_dJ )
                            {
                                bj_da_dJ = d_J_epart;

                                i2      = max_sec_pos[m];
                                e2      = max_sec_val[m];
                                oldtau2 = max_sec_old_tau[m];
                                newtau2 = max_sec_new_tau[m];

                                d_alpha1 = d_alpha1_trial;
                                d_alpha2 = d_alpha2_trial;
                                d_b      = d_b_trial;
                                d_f      = d_f_trial;
                            }
                        }
                    }
                }

                if ( !first_try )
                {
                    step_taken = actually_take_step_d2c_pattern(bj_da_dJ,i1,i2,newtau1,newtau2,d_alpha1,d_alpha2,d_b,d_f,f_zero,f_zero_next,problem,sol_tol,J);
                }
            }
        }

        if ( !step_taken )
        {
            for ( m = 1 ; ( ( m <= problem.N ) && !step_taken ) ; m++ )
            {
                switch ( (problem.fast_tau)[m-1] )
                {
                    case 0:
                    {
                        switch ( (problem.fast_contype)[m-1] )
                        {
                            case 0: { break; }
                            case 1: { step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J); break; }
                            case 2: { step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J); break; }

                            default:
                            {
                                if ( rand() % 2 )
                                {
                                    if ( !( step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J) ) )
                                    {
                                        step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J);
                                    }
                                }

                                else
                                {
                                    if ( !( step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J) ) )
                                    {
                                        step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J);
                                    }
                                }

                                break;
                            }
                        }

                        break;
                    }

                    case -1: { step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J); break; }
                    case +1: { step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J); break; }
                    case +2: { step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J); break; }
                    default: { step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J); break; }
                }
            }
        }

        if ( !step_taken )
        {
            goto exit_point;
        }
    }
  }

  else
  {
    #ifdef DEBUG
    E_CHECKER_GEN
    #endif

    for ( ; ( ( (*loop_count) <= epochs ) && ( !(*async_exit_flag) ) ) ; (*loop_count)++ )
    {
        for ( i = 1 ; i <= problem.d2c_buff_size ; i++ )
        {
            max_pos_pos[i-1] = 0;
            max_neg_pos[i-1] = 0;

            max_pos_val[i-1] = 0.0;
            max_neg_val[i-1] = 0.0;
        }

        pos_buffer_len = 0;
        neg_buffer_len = 0;

        for ( i = 1 ; i <= problem.N ; i++ )
        {
            switch ( (problem.fast_tau)[i-1] )
            {
                case +2: { ee = problem.fast_e_tau[i-1]; SORTINTOPOS(ee,+2,+1); break; }
                case -2: { ee = problem.fast_e_tau[i-1]; SORTINTONEG(ee,-2,-1); break; }
                case +1: { ee = problem.fast_e_tau[i-1]; SORTINTOPOS(ee,+1,+1);
                                                         SORTINTONEG(ee,+1,+1); break; }
                case -1: { ee = problem.fast_e_tau[i-1]; SORTINTOPOS(ee,-1,-1);
                                                         SORTINTONEG(ee,-1,-1); break; }

                default:
                {
                    switch ( (problem.fast_contype)[i-1] )
                    {
                        case 0:  { break; }
                        case 1:  { ee = problem.fast_e_tau[i-1]; SORTINTOPOS(ee,0,-1); break; }
                        case 2:  { ee = problem.fast_e_tau[i-1]; SORTINTONEG(ee,0,+1); break; }
                        default: { ee = problem.fast_e_tau[i-1]; SORTINTOPOS(ee,0,-1);
                                                                 SORTINTONEG(ee,0,+1); break; }
                    }

                    break;
                }
            }
        }

        if ( f_zero && ( !(max_pos_pos[0]) || ( max_pos_val[0] <  sol_tol ) )
                    && ( !(max_neg_pos[0]) || ( max_neg_val[0] > -sol_tol ) ) )
        {
            goto exit_point;
        }

        step_taken = 0;

        i = 0;
        j = 0;

        while ( ( !step_taken        ) &&
                ( i < pos_buffer_len ) &&
                ( j < neg_buffer_len )    )
        {
            if ( max_pos_val[i] > -max_neg_val[j] )
            {
                i1      = max_pos_pos[i];
                e1      = max_pos_val[i];
                oldtau1 = max_pos_old_tau[i];
                newtau1 = max_pos_new_tau[i];

                i++;

                max_sec_pos        = max_neg_pos;
                max_sec_val        = max_neg_val;
                max_sec_new_tau    = max_neg_new_tau;
                max_sec_old_tau    = max_neg_new_tau;

                k = j;
                l = neg_buffer_len;
            }

            else
            {
                i1      = max_neg_pos[j];
                e1      = max_neg_val[j];
                oldtau1 = max_neg_old_tau[j];
                newtau1 = max_neg_new_tau[j];

                j++;

                max_sec_pos        = max_pos_pos;
                max_sec_val        = max_pos_val;
                max_sec_new_tau    = max_pos_new_tau;
                max_sec_old_tau    = max_pos_new_tau;

                k = i;
                l = pos_buffer_len;
            }

            if ( ( oldtau1 == -2 ) ||
                 ( oldtau1 ==  0 ) ||
                 ( oldtau1 == +2 )    )
            {
                for ( m = k ; ( ( m < l ) && !step_taken ) ; m++ )
                {
                    if ( ( max_sec_old_tau[m] == -2 ) ||
                         ( max_sec_old_tau[m] ==  0 ) ||
                         ( max_sec_old_tau[m] == +2 )    )
                    {
                        step_taken = takeStep_d2c_pattern_fastkern(i1,max_sec_pos[m],newtau1,max_sec_new_tau[m],e1,max_sec_val[m],f_zero,problem,sol_tol,J);
                    }
                }
            }

            if ( !step_taken )
            {
                first_try = 1;

                for ( m = k ; m < l ; m++ )
                {
                    if ( i1 != max_sec_pos[m] )
                    {
                        if ( trial_step_d2c_fastkern(d_J_epart,i1,max_sec_pos[m],newtau1,max_sec_new_tau[m],e1,max_sec_val[m],d_alpha1_trial,d_alpha2_trial,d_b_trial,d_f_trial,f_zero,f_zero_next,problem,sol_tol) )
                        {
                            if ( first_try )
                            {
                                bj_da_dJ = d_J_epart;

                                i2      = max_sec_pos[m];
                                e2      = max_sec_val[m];
                                oldtau2 = max_sec_old_tau[m];
                                newtau2 = max_sec_new_tau[m];

                                d_alpha1 = d_alpha1_trial;
                                d_alpha2 = d_alpha2_trial;
                                d_b      = d_b_trial;
                                d_f      = d_f_trial;

                                first_try = 0;
                            }

                            else if ( d_J_epart < bj_da_dJ )
                            {
                                bj_da_dJ = d_J_epart;

                                i2      = max_sec_pos[m];
                                e2      = max_sec_val[m];
                                oldtau2 = max_sec_old_tau[m];
                                newtau2 = max_sec_new_tau[m];

                                d_alpha1 = d_alpha1_trial;
                                d_alpha2 = d_alpha2_trial;
                                d_b      = d_b_trial;
                                d_f      = d_f_trial;
                            }
                        }
                    }
                }

                if ( !first_try )
                {
                    step_taken = actually_take_step_d2c_pattern(bj_da_dJ,i1,i2,newtau1,newtau2,d_alpha1,d_alpha2,d_b,d_f,f_zero,f_zero_next,problem,sol_tol,J);
                }
            }
        }

        if ( !step_taken )
        {
            for ( m = 1 ; ( ( m <= problem.N ) && !step_taken ) ; m++ )
            {
                switch ( (problem.fast_tau)[m-1] )
                {
                    case 0:
                    {
                        switch ( (problem.fast_contype)[m-1] )
                        {
                            case 0: { break; }
                            case 1: { step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J); break; }
                            case 2: { step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J); break; }

                            default:
                            {
                                if ( rand() % 2 )
                                {
                                    if ( !( step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J) ) )
                                    {
                                        step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J);
                                    }
                                }

                                else
                                {
                                    if ( !( step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J) ) )
                                    {
                                        step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J);
                                    }
                                }

                                break;
                            }
                        }

                        break;
                    }

                    case -1: { step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J); break; }
                    case +1: { step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J); break; }
                    case +2: { step_taken = examineExample_SMO(m,+1,f_zero,problem,sol_tol,J); break; }
                    default: { step_taken = examineExample_SMO(m,-1,f_zero,problem,sol_tol,J); break; }
                }
            }
        }

        if ( !step_taken )
        {
            goto exit_point;
        }
    }
  }

  exit_point:

  #ifdef DEBUG
  for ( i = 1 ; i <= problem.N ; i++ )
  {
    j = 0;

    switch ( problem.fast_tau[i-1] )
    {
        case -2: { if (      problem.e_tau[i-1]  >= 0.0  ) { j = 1; } break; }
        case +2: { if (      problem.e_tau[i-1]  <= 0.0  ) { j = 1; } break; }
        case -1: { if ( fabs(problem.e_tau[i-1]) <= 1e06 ) { j = 1; } break; }
        case +1: { if ( fabs(problem.e_tau[i-1]) <= 1e06 ) { j = 1; } break; }
        default:
        {
            if ( problem.fast_contype[i-1] == 2 )
            {
                if ( problem.e_tau[i-1] >= 0.0 )
                {
                    j = 1;
                }
            }

            else
            {
                if ( problem.e_tau[i-1] <= 0.0 )
                {
                    j = 1;
                }
            }
        }
    }

    if ( j == 0 )
    {
      std::cerr << "point " << i << " is not optimal\n";
    }
  }
  #endif

  REPDELB(max_pos_pos);
  REPDELB(max_pos_val);
  REPDELB(max_pos_old_tau);
  REPDELB(max_pos_new_tau);
  REPDELB(max_neg_pos);
  REPDELB(max_neg_val);
  REPDELB(max_neg_old_tau);
  REPDELB(max_neg_new_tau);


  return (*loop_count);
}

unsigned long solve_d2c_fixed_bias(SVdata &problem, L_DOUBLE &sol_tol, unsigned long epochs, volatile int *async_exit_flag, volatile unsigned long *loop_count)
{
    /*
       D2C undefined for this, but logically would reduce to SMO (which also
       doesn't exist anywhere but here, as an obvious extension of Platt's
       original variable bias algorithm.
    */

    return solve_SMO_fixed_bias(problem,sol_tol,epochs,async_exit_flag,loop_count);
}












int takeStep_d2c(long i1, long i2, int tau1, int tau2, L_DOUBLE &e1, L_DOUBLE &e2, int f_zero, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J)
{
    if ( i1 == i2 )
    {
        return 0;
    }

    L_DOUBLE d_alpha1 = 0.0;
    L_DOUBLE d_alpha2 = 0.0;
    L_DOUBLE d_b = 0.0;
    L_DOUBLE d_f = 0.0;
    L_DOUBLE d_J_epart = 0.0;
    int f_zero_next = f_zero;

    if ( trial_step_d2c(d_J_epart,i1,i2,tau1,tau2,e1,e2,d_alpha1,d_alpha2,d_b,d_f,f_zero,f_zero_next,problem,sol_tol) == 0 )
    {
        return 0;
    }

    return actually_take_step_d2c(d_J_epart,i1,i2,tau1,tau2,d_alpha1,d_alpha2,d_b,d_f,f_zero,f_zero_next,problem,sol_tol,J);
}

int takeStep_d2c_pattern(long i1, long i2, int tau1, int tau2, L_DOUBLE &e1, L_DOUBLE &e2, int f_zero, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J)
{
    if ( i1 == i2 )
    {
        return 0;
    }

    L_DOUBLE d_alpha1 = 0.0;
    L_DOUBLE d_alpha2 = 0.0;
    L_DOUBLE d_b = 0.0;
    L_DOUBLE d_f = 0.0;
    L_DOUBLE d_J_epart = 0.0;
    int f_zero_next = f_zero;

    if ( trial_step_d2c(d_J_epart,i1,i2,tau1,tau2,e1,e2,d_alpha1,d_alpha2,d_b,d_f,f_zero,f_zero_next,problem,sol_tol) == 0 )
    {
        return 0;
    }

    return actually_take_step_d2c_pattern(d_J_epart,i1,i2,tau1,tau2,d_alpha1,d_alpha2,d_b,d_f,f_zero,f_zero_next,problem,sol_tol,J);
}

int takeStep_d2c_pattern_fastkern(long i1, long i2, int tau1, int tau2, L_DOUBLE &e1, L_DOUBLE &e2, int f_zero, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J)
{
    if ( i1 == i2 )
    {
        return 0;
    }

    L_DOUBLE d_alpha1 = 0.0;
    L_DOUBLE d_alpha2 = 0.0;
    L_DOUBLE d_b = 0.0;
    L_DOUBLE d_f = 0.0;
    L_DOUBLE d_J_epart = 0.0;
    int f_zero_next = f_zero;

    if ( trial_step_d2c_fastkern(d_J_epart,i1,i2,tau1,tau2,e1,e2,d_alpha1,d_alpha2,d_b,d_f,f_zero,f_zero_next,problem,sol_tol) == 0 )
    {
        return 0;
    }

    return actually_take_step_d2c_pattern(d_J_epart,i1,i2,tau1,tau2,d_alpha1,d_alpha2,d_b,d_f,f_zero,f_zero_next,problem,sol_tol,J);
}

int trial_step_d2c(L_DOUBLE &d_J_epart, long i1, long i2, int tau1, int tau2, L_DOUBLE &e1, L_DOUBLE &e2, L_DOUBLE &d_alpha1, L_DOUBLE &d_alpha2, L_DOUBLE &d_b, L_DOUBLE &d_f, int &f_zero, int &f_zero_next, SVdata &problem, L_DOUBLE &sol_tol)
{
    //
    // Now:
    //
    // [ d_b      ]          [ 0  1   1  ]   [  f  ]
    // [ d_alpha1 ] = - inv( [ 1 K11 K12 ] ) [ e_1 ]
    // [ d_alpha2 ]          [ 1 K12 K22 ]   [ e_2 ]
    //
    //                       1            [  K11.K22 - K12.K12 ]
    //              = ----------------- ( [      K12 - K22     ] f
    //                K11 + K22 - 2.K12   [      K12 - K11     ]
    //
    //                  [ K12 - K22 ]       [ K12 - K11 ]
    //                + [    -1     ] e_1 + [    +1     ] e_2 )
    //                  [    +1     ]       [    -1     ]
    //
    // [ d_f  ]   [ -d_f  ]
    // [ d_e1 ] = [ -d_e1 ]
    // [ d_e2 ]   [ -d_e2 ]
    //
    // Of course, this is only true if the 2 hessian is non-singular - ie.
    // K11 + K22 - 2.K12 != 0.  If this does not hold, no problem - just
    // find a direction of linear non-ascent wrt alpha_1 and alpha_2,
    // ignoring f in this case.
    //
    // In the singular case:
    //
    // [ d_b      ] = - theta inv( [ 0  1  ] ) [  1  ]
    // [ d_alpha1 ]                [ 1 K11 ]   [ K12 ]
    //
    //              = - theta (1/-1) [ K11 -1 ] [  1  ]
    //                               [ -1  0  ] [ K12 ]
    //
    //              = theta [ K11 -1 ] [  1  ]
    //                      [ -1  0  ] [ K12 ]
    //
    //              = theta [ K11 - K12 ]
    //                      [    -1     ]
    //
    // d_alpha2 = theta
    //
    // [ d_f  ]   [ 0 ]
    // [ d_e1 ] = [ 0 ]
    // [ d_e2 ]   [ 0 ]
    //
    // where theta is chosen to ensure linear non-increase (e1.d_alpha1 +
    // e2.d_alpha2 < 0) and also to make sure that a bound is hit.
    // Specifically:
    //
    // theta = -1.1 * |h1| : if e2 >= e1 and alpha1 is positive
    // theta = -1.1 * |v1| : if e2 >= e1 and alpha1 is negative
    // theta = +1.1 * |h1| : if e2 <  e1 and alpha1 is positive
    // theta = +1.1 * |v1| : if e2 <  e1 and alpha1 is negative
    //

    L_DOUBLE K11;
    L_DOUBLE K22;
    L_DOUBLE K12;
    L_DOUBLE negdetH;
    L_DOUBLE f;

    K11 = problem.opt_get_G_tau(i1-1,i1-1);
    K22 = problem.opt_get_G_tau(i2-1,i2-1);
    K12 = problem.opt_get_G_tau(i1-1,i2-1);

    negdetH = K11+K22-(2.0*K12);

    f = problem.f;

    if ( ( negdetH > ZZERO ) || ( negdetH < -ZZERO ) )
    {
        /*
           Nonsingular case
        */

        if ( f_zero )
        {
            d_b = (((K12-K22)*e1)+((K12-K11)*e2))/negdetH;
            d_f = 0.0;

            d_alpha1 = (-e1+e2)/negdetH;
            d_alpha2 = (+e1-e2)/negdetH;

            d_J_epart = -((e1-e2)*(e1-e2))/negdetH;
        }

        else
        {
            d_b = ((((K11*K22)-(K12*K12))*f)+((K12-K22)*e1)+((K12-K11)*e2))/negdetH;
            d_f = -(problem.f);

            d_alpha1 = (((K12-K22)*f)-e1+e2)/negdetH;
            d_alpha2 = (((K12-K11)*f)+e1-e2)/negdetH;

            d_J_epart = -((e1-e2)*(e1-e2))/negdetH;
        }
    }

    else
    {
        /*
           Singular case
        */

        L_DOUBLE theta;

        d_b = K11-K12;
        d_f = 0.0;

        d_alpha1 = -1.0;
        d_alpha2 = +1.0;

        if ( tau2 == +1 )
        {
            if ( e2 >= e1 )
            {
                theta = -1.1 * (problem.fast_h)[i2];
            }

            else
            {
                theta = 1.1 * (problem.fast_h)[i2];
            }
        }

        else
        {
            if ( e2 >= e1 )
            {
                theta = -1.1 * -(problem.fast_v)[i2];
            }

            else
            {
                theta = 1.1 * -(problem.fast_v)[i2];
            }
        }

        d_b *= theta;

        d_alpha1 *= theta;
        d_alpha2 *= theta;

        d_J_epart = (d_alpha1*e1) + (d_alpha2*e2);
    }

    /*
       Now let's scale the step
    */

    L_DOUBLE step_scale = 1.0;
    L_DOUBLE ceta;

    f_zero_next = 1;

    if ( ( d_alpha1 > -ZZERO/10 ) && ( d_alpha1 < ZZERO/10 ) )
    {
        d_alpha1 = 0.0;
    }

    if ( ( d_alpha2 > -ZZERO/10 ) && ( d_alpha2 < ZZERO/10 ) )
    {
        d_alpha2 = 0.0;
    }

    if ( d_alpha1 < 0.0 )
    {
        if ( tau1 > 0 )
        {
            ceta = -((problem.fast_alpha)[i1-1]) / d_alpha1;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }

        else
        {
            ceta = ( (problem.fast_v)[i1-1] - ((problem.fast_alpha)[i1-1]) ) / d_alpha1;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }
    }

    else if ( d_alpha1 > 0.0 )
    {
        if ( tau1 > 0 )
        {
            ceta = ( (problem.fast_h)[i1-1] - ((problem.fast_alpha)[i1-1]) ) / d_alpha1;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }

        else
        {
            ceta = -((problem.fast_alpha)[i1-1]) / d_alpha1;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }
    }

    if ( d_alpha2 < 0.0 )
    {
        if ( tau2 > 0 )
        {
            ceta = -((problem.fast_alpha)[i2-1]) / d_alpha2;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }

        else
        {
            ceta = ( (problem.fast_v)[i2-1] - ((problem.fast_alpha)[i2-1]) ) / d_alpha2;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }
    }

    else if ( d_alpha2 > 0.0 )
    {
        if ( tau2 > 0 )
        {
            ceta = ( (problem.fast_h)[i2-1] - ((problem.fast_alpha)[i2-1]) ) / d_alpha2;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }

        else
        {
            ceta = -((problem.fast_alpha)[i2-1]) / d_alpha2;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }
    }

    d_b *= step_scale;
    d_f *= step_scale;

    d_alpha1 *= step_scale;
    d_alpha2 *= step_scale;

    if ( ( negdetH > ZZERO ) || ( negdetH < -ZZERO ) )
    {
        /*
           Step in quadratically curved direction, so scaling...
        */

        d_J_epart *= (step_scale-((step_scale*step_scale)/2.0));
    }

    else
    {
        /*
           Step in linear direction, so scaling is linear.
        */

        d_J_epart *= step_scale;
    }

    /*
       Return zero if step too small.  This is from Platt, but can only
       be applied validly if f is zero (as it assumes that d_alpha1 and
       d_alpha2 are the same up to sign).
    */

    if ( f_zero && ( fabs(d_alpha1) < sol_tol*(fabs((2.0*((problem.fast_alpha)[i1-1]))+d_alpha1)+sol_tol) ) )
    {
        return 0;
    }

    return 1;
}

int trial_step_d2c_fastkern(L_DOUBLE &d_J_epart, long i1, long i2, int tau1, int tau2, L_DOUBLE &e1, L_DOUBLE &e2, L_DOUBLE &d_alpha1, L_DOUBLE &d_alpha2, L_DOUBLE &d_b, L_DOUBLE &d_f, int &f_zero, int &f_zero_next, SVdata &problem, L_DOUBLE &sol_tol)
{
    //
    // Now:
    //
    // [ d_b      ]          [ 0  1   1  ]   [  f  ]
    // [ d_alpha1 ] = - inv( [ 1 K11 K12 ] ) [ e_1 ]
    // [ d_alpha2 ]          [ 1 K12 K22 ]   [ e_2 ]
    //
    //                       1            [  K11.K22 - K12.K12 ]
    //              = ----------------- ( [      K12 - K22     ] f
    //                K11 + K22 - 2.K12   [      K12 - K11     ]
    //
    //                  [ K12 - K22 ]       [ K12 - K11 ]
    //                + [    -1     ] e_1 + [    +1     ] e_2 )
    //                  [    +1     ]       [    -1     ]
    //
    // [ d_f  ]   [ -d_f  ]
    // [ d_e1 ] = [ -d_e1 ]
    // [ d_e2 ]   [ -d_e2 ]
    //
    // Of course, this is only true if the 2 hessian is non-singular - ie.
    // K11 + K22 - 2.K12 != 0.  If this does not hold, no problem - just
    // find a direction of linear non-ascent wrt alpha_1 and alpha_2,
    // ignoring f in this case.
    //
    // In the singular case:
    //
    // [ d_b      ] = - theta inv( [ 0  1  ] ) [  1  ]
    // [ d_alpha1 ]                [ 1 K11 ]   [ K12 ]
    //
    //              = - theta (1/-1) [ K11 -1 ] [  1  ]
    //                               [ -1  0  ] [ K12 ]
    //
    //              = theta [ K11 -1 ] [  1  ]
    //                      [ -1  0  ] [ K12 ]
    //
    //              = theta [ K11 - K12 ]
    //                      [    -1     ]
    //
    // d_alpha2 = theta
    //
    // [ d_f  ]   [ 0 ]
    // [ d_e1 ] = [ 0 ]
    // [ d_e2 ]   [ 0 ]
    //
    // where theta is chosen to ensure linear non-increase (e1.d_alpha1 +
    // e2.d_alpha2 < 0) and also to make sure that a bound is hit.
    // Specifically:
    //
    // theta = -1.1 * |h1| : if e2 >= e1 and alpha1 is positive
    // theta = -1.1 * |v1| : if e2 >= e1 and alpha1 is negative
    // theta = +1.1 * |h1| : if e2 <  e1 and alpha1 is positive
    // theta = +1.1 * |v1| : if e2 <  e1 and alpha1 is negative
    //

    L_DOUBLE K11;
    L_DOUBLE K22;
    L_DOUBLE K12;
    L_DOUBLE negdetH;
    L_DOUBLE f;

    K11 = (problem.fast_kernel_fn)((problem.fast_opt_x)[i1-1],(problem.fast_opt_x)[i1-1],problem.fast_uc,problem.fast_uv,(double *) problem.fast_covw);
    K22 = (problem.fast_kernel_fn)((problem.fast_opt_x)[i2-1],(problem.fast_opt_x)[i2-1],problem.fast_uc,problem.fast_uv,(double *) problem.fast_covw);
    K12 = (problem.fast_kernel_fn)((problem.fast_opt_x)[i1-1],(problem.fast_opt_x)[i2-1],problem.fast_uc,problem.fast_uv,(double *) problem.fast_covw);

    negdetH = K11+K22-(2.0*K12);

    f = problem.f;

    if ( ( negdetH > ZZERO ) || ( negdetH < -ZZERO ) )
    {
        /*
           Nonsingular case
        */

        if ( f_zero )
        {
            d_b = (((K12-K22)*e1)+((K12-K11)*e2))/negdetH;
            d_f = 0.0;

            d_alpha1 = (-e1+e2)/negdetH;
            d_alpha2 = (+e1-e2)/negdetH;

            d_J_epart = -((e1-e2)*(e1-e2))/negdetH;
        }

        else
        {
            d_b = ((((K11*K22)-(K12*K12))*f)+((K12-K22)*e1)+((K12-K11)*e2))/negdetH;
            d_f = -(problem.f);

            d_alpha1 = (((K12-K22)*f)-e1+e2)/negdetH;
            d_alpha2 = (((K12-K11)*f)+e1-e2)/negdetH;

            d_J_epart = -((e1-e2)*(e1-e2))/negdetH;
        }
    }

    else
    {
        /*
           Singular case
        */

        L_DOUBLE theta;

        d_b = K11-K12;
        d_f = 0.0;

        d_alpha1 = -1.0;
        d_alpha2 = +1.0;

        if ( tau2 == +1 )
        {
            if ( e2 >= e1 )
            {
                theta = -1.1 * (problem.fast_h)[i2];
            }

            else
            {
                theta = 1.1 * (problem.fast_h)[i2];
            }
        }

        else
        {
            if ( e2 >= e1 )
            {
                theta = -1.1 * -(problem.fast_v)[i2];
            }

            else
            {
                theta = 1.1 * -(problem.fast_v)[i2];
            }
        }

        d_b *= theta;

        d_alpha1 *= theta;
        d_alpha2 *= theta;

        d_J_epart = (d_alpha1*e1) + (d_alpha2*e2);
    }

    /*
       Now let's scale the step
    */

    L_DOUBLE step_scale = 1.0;
    L_DOUBLE ceta;

    f_zero_next = 1;

    if ( ( d_alpha1 > -ZZERO/10 ) && ( d_alpha1 < ZZERO/10 ) )
    {
        d_alpha1 = 0.0;
    }

    if ( ( d_alpha2 > -ZZERO/10 ) && ( d_alpha2 < ZZERO/10 ) )
    {
        d_alpha2 = 0.0;
    }

    if ( d_alpha1 < 0.0 )
    {
        if ( tau1 > 0 )
        {
            ceta = -((problem.fast_alpha)[i1-1]) / d_alpha1;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }

        else
        {
            ceta = ( (problem.fast_v)[i1-1] - ((problem.fast_alpha)[i1-1]) ) / d_alpha1;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }
    }

    else if ( d_alpha1 > 0.0 )
    {
        if ( tau1 > 0 )
        {
            ceta = ( (problem.fast_h)[i1-1] - ((problem.fast_alpha)[i1-1]) ) / d_alpha1;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }

        else
        {
            ceta = -((problem.fast_alpha)[i1-1]) / d_alpha1;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }
    }

    if ( d_alpha2 < 0.0 )
    {
        if ( tau2 > 0 )
        {
            ceta = -((problem.fast_alpha)[i2-1]) / d_alpha2;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }

        else
        {
            ceta = ( (problem.fast_v)[i2-1] - ((problem.fast_alpha)[i2-1]) ) / d_alpha2;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }
    }

    else if ( d_alpha2 > 0.0 )
    {
        if ( tau2 > 0 )
        {
            ceta = ( (problem.fast_h)[i2-1] - ((problem.fast_alpha)[i2-1]) ) / d_alpha2;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }

        else
        {
            ceta = -((problem.fast_alpha)[i2-1]) / d_alpha2;

            if ( ceta < step_scale )
            {
                step_scale = ceta;
                f_zero_next = f_zero;
            }
        }
    }

    d_b *= step_scale;
    d_f *= step_scale;

    d_alpha1 *= step_scale;
    d_alpha2 *= step_scale;

    if ( ( negdetH > ZZERO ) || ( negdetH < -ZZERO ) )
    {
        /*
           Step in quadratically curved direction, so scaling...
        */

        d_J_epart *= (step_scale-((step_scale*step_scale)/2.0));
    }

    else
    {
        /*
           Step in linear direction, so scaling is linear.
        */

        d_J_epart *= step_scale;
    }

    /*
       Return zero if step too small.  This is from Platt, but can only
       be applied validly if f is zero (as it assumes that d_alpha1 and
       d_alpha2 are the same up to sign).
    */

    if ( f_zero && ( fabs(d_alpha1) < sol_tol*(fabs((2.0*((problem.fast_alpha)[i1-1]))+d_alpha1)+sol_tol) ) )
    {
        return 0;
    }

    return 1;
}

int actually_take_step_d2c(L_DOUBLE &d_J_epart, long i1, long i2, int tau1, int tau2, L_DOUBLE &d_alpha1, L_DOUBLE &d_alpha2, L_DOUBLE &d_b, L_DOUBLE &d_f, int &f_zero, int &f_zero_next, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J)
{
    /*
       ...but first, we may need to unconstrain alpha1 and/or alpha2
    */

    if ( ( (problem.fast_tau)[i1-1] != +1 ) && ( (problem.fast_tau)[i1-1] != -1 ) )
    {
        if ( tau1 == +1 )
        {
            (problem.fast_tau)[i1-1] = +1;
        }

        else
        {
            (problem.fast_tau)[i1-1] = -1;
        }
    }

    if ( ( (problem.fast_tau)[i2-1] != +1 ) && ( (problem.fast_tau)[i2-1] != -1 ) )
    {
        if ( tau2 == +1 )
        {
            (problem.fast_tau)[i2-1] = +1;
        }

        else
        {
            (problem.fast_tau)[i2-1] = -1;
        }
    }

    (*J) += d_J_epart;

    problem.opt_step_two(d_b,d_f,d_alpha1,d_alpha2,i1-1,i2-1);

    f_zero = f_zero_next;

    /*
       Be careful with constraints.  If alpha wanders too close to the
       boundary, reel it in and constrain it.
    */

    if ( ( (problem.fast_alpha)[i1-1] <= ZZERO ) && ( (problem.fast_alpha)[i1-1] >= -ZZERO ) )
    {
        (problem.fast_alpha)[i1-1] = 0.0;

        problem.opt_constrain_Z(i1-1);
    }

    else if ( ( (problem.fast_alpha)[i1-1] <= ((problem.fast_v)[i1-1])+ZZERO ) && ( (problem.fast_tau)[i1-1] == -1 ) )
    {
        (problem.fast_alpha)[i1-1] = (problem.fast_v)[i1-1];
        (problem.fast_tau)[i1-1]   = -2;
    }

    else if ( ( (problem.fast_alpha)[i1-1] >= ((problem.fast_h)[i1-1])-ZZERO ) && ( (problem.fast_tau)[i1-1] == +1 ) )
    {
        (problem.fast_alpha)[i1-1] = (problem.fast_h)[i1-1];
        (problem.fast_tau)[i1-1]   = +2;
    }

    if ( ( (problem.fast_alpha)[i2-1] <= ZZERO ) && ( (problem.fast_alpha)[i2-1] >= -ZZERO ) )
    {
        (problem.fast_alpha)[i2-1] = 0.0;

        problem.opt_constrain_Z(i2-1);
    }

    else if ( ( (problem.fast_alpha)[i2-1] <= ((problem.fast_v)[i2-1])+ZZERO ) && ( (problem.fast_tau)[i2-1] == -1 ) )
    {
        (problem.fast_alpha)[i2-1] = (problem.fast_v)[i2-1];
        (problem.fast_tau)[i2-1]   = -2;
    }

    else if ( ( (problem.fast_alpha)[i2-1] >= ((problem.fast_h)[i2-1])-ZZERO ) && ( (problem.fast_tau)[i2-1] == +1 ) )
    {
        (problem.fast_alpha)[i2-1] = (problem.fast_h)[i2-1];
        (problem.fast_tau)[i2-1]   = +2;
    }

    return 1;

    /* burn warnings burn */

    sol_tol = 0.0;
}

int actually_take_step_d2c_pattern(L_DOUBLE &d_J_epart, long i1, long i2, int tau1, int tau2, L_DOUBLE &d_alpha1, L_DOUBLE &d_alpha2, L_DOUBLE &d_b, L_DOUBLE &d_f, int &f_zero, int &f_zero_next, SVdata &problem, L_DOUBLE &sol_tol, JTYPE *J)
{
    /*
       ...but first, we may need to unconstrain alpha1 and/or alpha2
    */

    if ( ( (problem.fast_tau)[i1-1] != +1 ) && ( (problem.fast_tau)[i1-1] != -1 ) )
    {
        if ( tau1 == +1 )
        {
            (problem.fast_tau)[i1-1] = +1;
        }

        else
        {
            (problem.fast_tau)[i1-1] = -1;
        }
    }

    if ( ( (problem.fast_tau)[i2-1] != +1 ) && ( (problem.fast_tau)[i2-1] != -1 ) )
    {
        if ( tau2 == +1 )
        {
            (problem.fast_tau)[i2-1] = +1;
        }

        else
        {
            (problem.fast_tau)[i2-1] = -1;
        }
    }

    (*J) += d_J_epart;

    problem.opt_step_two(d_b,d_f,d_alpha1,d_alpha2,i1-1,i2-1);

    f_zero = f_zero_next;

    /*
       Be careful with constraints.  If alpha wanders too close to the
       boundary, reel it in and constrain it.
    */

    if ( ( (problem.fast_alpha)[i1-1] <= ZZERO ) && ( (problem.fast_alpha)[i1-1] >= -ZZERO ) )
    {
        (problem.fast_alpha)[i1-1] = 0.0;
        (problem.fast_tau)[i1-1]   = 0;
    }

    else if ( ( (problem.fast_alpha)[i1-1] <= ((problem.fast_v)[i1-1])+ZZERO ) && ( (problem.fast_tau)[i1-1] == -1 ) )
    {
        (problem.fast_alpha)[i1-1] = (problem.fast_v)[i1-1];
        (problem.fast_tau)[i1-1]   = -2;
    }

    else if ( ( (problem.fast_alpha)[i1-1] >= ((problem.fast_h)[i1-1])-ZZERO ) && ( (problem.fast_tau)[i1-1] == +1 ) )
    {
        (problem.fast_alpha)[i1-1] = (problem.fast_h)[i1-1];
        (problem.fast_tau)[i1-1]   = +2;
    }

    if ( ( (problem.fast_alpha)[i2-1] <= ZZERO ) && ( (problem.fast_alpha)[i2-1] >= -ZZERO ) )
    {
        (problem.fast_alpha)[i2-1] = 0.0;
        (problem.fast_tau)[i2-1]   = 0;
    }

    else if ( ( (problem.fast_alpha)[i2-1] <= ((problem.fast_v)[i2-1])+ZZERO ) && ( (problem.fast_tau)[i2-1] == -1 ) )
    {
        (problem.fast_alpha)[i2-1] = (problem.fast_v)[i2-1];
        (problem.fast_tau)[i2-1]   = -2;
    }

    else if ( ( (problem.fast_alpha)[i2-1] >= ((problem.fast_h)[i2-1])-ZZERO ) && ( (problem.fast_tau)[i2-1] == +1 ) )
    {
        (problem.fast_alpha)[i2-1] = (problem.fast_h)[i2-1];
        (problem.fast_tau)[i2-1]   = +2;
    }

    return 1;

    /* burn warnings burn */

    sol_tol = 0.0;
}











