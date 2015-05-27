#include "stdafx.h"
/*
 *  RIMElib: RuntIme Mathematical Equation Library
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

#include <math.h>
#include <float.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#include <gsl/gsl_math.h>
#include <gsl/specfunc/gsl_sf.h>
#include <gsl/err/gsl_errno.h>

#include "gvars.h"
#include "gmaths.h"
#include "erime.h"
#include "mathtext.h"

#ifndef INIT_EMATHS
#include "erime.dat"
#endif

/*********************************************************************

Basic data types
================

MathsNode: Equations are stored as a tree, where each vertex of this
        tree is of type MathsNode.  MathsNode contains the following:

        CombinerFunction: the function number that this represents.
        multiplier:     scale factor for result.
        SubFunctions:   array of pointers to arguments used by this
                        function (ie. pointers to the next layer
                        down the tree).

GenArg: What gets passed about when evaluating a function.  What is
        stored here can be divided into local and global variables.

        The global variables are a (void) pointer to some generic
        "data" (argContents) and a function that can make sense of
        this argument (ArgEvaluationFunction).  When var(i,j), i>0,
        j>0 is evaluated, this function will be called and passed
        argContents, i and j.  It (should) return either some value
        (possibly an error if i or j are out of range) or
        DATA_IS_UNCALCED if data is not available.

        The following functions are the "standard" functions for
        getting values of variables:

        ArgEvalNull:    returns DATA_IS_UNCALCED.  This is used when
                        simplifying functions to indicate that this
                        cannot be simplified out (whereas a function
                        that returns a definite number may be).
        ArgEvalDouble:  return element in 1d array (double *)
                        argContents where j is the index and i must be
                        1.  Otherwise, it returns an error.
        ArgEvalDoubleDouble: return element in 2d array (double **)
                        argContents.
        ArgEvalList:    like ArgEvalDouble, but with bounds checking
                        to ensure that i is not too large.  If oob
                        will return an error.
        ArgEvalGeneral: given (GenArgGeneral *) argContents will call
                        another function (which is external to this
                        library) to get the result.

        "Data" for these functions may be stored in the following
        structures:

        GenArgList:     contains a pointer to an array and the size of
                        that array.
        GenArgGeneral:  contains a pointer to some data and a pointer
                        to a function capable of deciphering it, given
                        that pointer and also the indices i and j.
                        This is just like GenArg, but usable
                        externally.

        Local variables are allocated when needed and contain "local"
        variables - ie. those refered to as var(0,i), i>0, which may
        appear in the equation as summation indices etc.  The data
        is kept in loc_vars[i-1] (size MAX_NUM_LOC_VARS).  
        loc_var_used_flags[i-1] is 1 if var(0,i) appears (or possibly
        has appeared but is no longer present) in some equation and 0
        otherwise.  This must be allocated before use (see below).

FunctionNode: General function descriptor struct - this holds details
        about a particular function, namely:

        FnName:         a string with the name of this function, for
                        parsing purposes.
        NumLocInput:    number of arguments taken by this function.
        CalcOperation:  function to evaluate this function.
        Deriv_string:   this string holds the derivative of the
                        function.  It is parsed by initMaths to create
                        static data Deriv_equ.  var(i,j) where i and j
                        are negative are special variables that must
                        be replaced when constructing the the actual
                        derivative for a concrete case.  There are:
                        var(-1,-i): SubFunction (argument) i
                        var(-2,-i): Derivative of SubFunction (arg) i
                        var(-3,-i): Local variable placeholder.  This
                                will be replaced with a number that
                                represents a local variable that is
                                not used by either this node or its
                                subfunctions.  Here, 1 <= i <=
                                NumLocDerivVar
                        var(-4,-1): This will be replaced with i where
                                deriv is w.r.t. var(i,j).
                        var(-4,-2): This will be replaced with j where
                                deriv is w.r.t. var(i,j).
                        Note:   for reasons of code neatness, these
                                derivatives are declared in a bunch of
                                strings of the form deriv_***.
                                Pointers to these are then included in
                                this structure.
        Deriv_equ:      Parsed form of Deriv_string. Stored statically
                        in erime.dat.  This file is constructed using
                        the function initMaths.
        NumLocDerivVar: Number of local variable placeholders that
                        need to be replaced.
        AltMarker:      usually ',', but if another value is used then
                        during deparsing FnName will not be appended
                        to the start of this function.  Rather, the
                        separator between arguments (which is normally
                        ,) will be replaced by this.  For example, the
                        function add uses AltMarker = '+' to create
                        human readable output.

*********************************************************************/

typedef struct
{
    void *argContents;

    int (*ArgEvaluationFunction)(GenVar *result, void *argContents, long i, long j);

    GenVar **loc_vars;
    int **loc_var_used_flags;
}
GenArg;

typedef struct
{
    int size;

    double *argContents;
}
GenArgList;

typedef struct
{
    void *argContents;

    double (*ArgEvaluationFunction)(void *argContents, long i, long j);
}
GenArgGeneral;

typedef struct
{
    char *FnName;
    int NumLocInput;

    int (*CalcOperation)(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);

    char *Deriv_string;
    MathsNode *Deriv_equ;
    int NumLocDerivVar;

    char altMarker;
}
FunctionNode;

int ArgEvalNull(GenVar *result, void *argContents, long i, long j);
int ArgEvalDouble(GenVar *result, void *argContents, long i, long j);
int ArgEvalDoubleDouble(GenVar *result, void *argContents, long i, long j);
int ArgEvalList(GenVar *result, void *argContents, long i, long j);
int ArgEvalGeneral(GenVar *result, void *argContents, long i, long j);






/*********************************************************************

Internal Functions
==================

The functions prefixed by l (on its own) are internal versions of the
external functions, differing only in the lack of post-simplification
(so the result of lparseMaths, for example, is exactly what you put
in, even if that evaluates to a constant).

Other functions are:

locparseMaths:  parse equation assumed to be in most basic form (ie.
                all brackets standardised, all operators removed etc.)
locdeparseMaths_cstruct: Subfunction used by deparseMaths_cstruct.

locsubVar:      just like fast_subVar, but without bounds checking.
locsubVarRow:   just like fast_subVarRow, but without bounds checking.
locsubVarCol:   just like fast_subVarCol, but without bounds checking.

locmakeDeriv:   make the derivative.  The important thing about this
                is the presence of GlobInput, which is necessary as
                the local variables are global in the context of this
                equation, so we must keep track of what local varialbes
                are available throughout this operation.
locMatrixMakeDeriv: The above function in the general case.
locRowMakeDeriv: The above function in the general case.
locColMakeDeriv: The above function in the general case.

locevaluateEqn: evaluate equation given argument.  The areguments are:
                touch_all - if non-zero, all branches of the tree will
                        be followed even if they don't strictly effect
                        the result.  This is needed when we want to
                        track what local variables will be available
                        for use - eg. while calculating a derivative.
                report_change - this affects the random number
                        generators only.  Rather than returning a
                        random number, if this is set they will return
                        "uncalced".  This is important in several
                        cases.  eg. simplification - if this is not
                        set, the simplify will simply convert the
                        random number to a constant, which is not what
                        is wanted.
                simplify_onthefly - Modifies the equation as evaluation
                        proceeds to make it smaller.  In particular,
                        if a node evaluates to a number (ie. known) it
                        will replace the node with this number.  If
                        this is set then, although local vars (and
                        flags) are set as usual, they will not be used.
                        That is, var(0,i) will evaluate to uncalculated.
                        This prevents the simplifier removing local
                        variables from the equation.
make_number:    If a function evaluates to a constant, then this
                constant (prior to multiplication by the multiplier)
                and a pointer to the function itself may be passed to
                this function.  This function will delete the equation
                (the subfunctions) and replace it with an appropriate
                constant.  It will then return an appropriate constant
                evaluation of this function such that when this is
                multiplied by the new multiplier it will give the same
                result as the function itself will (when the multiplier
                is multiplied by... my god, I think I'm turning into
                Dr Seus).

*********************************************************************/

MathsNode *locparseMaths(char *express, long len);
char *locdeparseMaths_cstruct(MathsNode *equation, char *name, int namelen, int *index);

int locsubVar(MathsNode *equation, MathsNode *var_replace, long i, long j, int naive);
int locsubVarRow(MathsNode *equation, MathsNode *var_replace, long i, int naive);
int locsubVarCol(MathsNode *equation, MathsNode *var_replace, long j, int naive);

MathsNode *locmakeDeriv(MathsNode *start, long i, long j, GenArg GlobInput);
MathsNode *locMatrixMakeDeriv(MathsNode *start, long i, long j, long k, long l, GenArg GlobInput);
MathsNode *locRowMakeDeriv(MathsNode *start, long i, long j, long k, GenArg GlobInput);
MathsNode *locColMakeDeriv(MathsNode *start, long i, long j, long k, GenArg GlobInput);

int locevaluateEqn(GenVar *result, MathsNode *equation, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);
int Gen_evaluateEqnNull(GenVar *result, MathsNode *equation, int report_change, int simplify_onthefly);
int Gen_evaluateEqnMatrix(GenVar *result, MathsNode *equation, double **GlobInput);
int Gen_evaluateEqnVector(GenVar *result, MathsNode *equation, double *GlobInput);
int Gen_evaluateEqnGeneral(GenVar *result, MathsNode *equation, double (*DoubleArgEvaluationFunction)(void *argContents, long i, long j), void *GlobInput);
void make_number(MathsNode *equation, GenVar *result);





/*********************************************************************

Derivative strings
==================

The following strings define the result of differentiating a function
of a particular type.  See above for details on contents.  Changing
these does not automatically change the derivatives - they are in the
erime.dat file.  To change a derivative, change the relevant string,
run initMaths to rebuild the static variables and cut/paste the result
into erime.dat (but make a backup first).

Various globals are defined - see initMaths for details.

*********************************************************************/

char Deriv_one_string[]               = "0";
char Deriv_zero_string[]              = "0";
char Deriv_inf_string[]               = "0";
char Deriv_ninf_string[]              = "0";
char Deriv_FiniteIndet_string[]       = "0";
char Deriv_InfintIndet_string[]       = "0";
char Deriv_isError_string[]           = "isError()";
char Deriv_var_string[]               = "kronDelta(var(-1,-1),var(-4,-1))*kronDelta(var(-1,-2),var(-4,-2))";
char Deriv_x_string[]                 = "kronDelta(1,var(-4,-1))*kronDelta(1,var(-4,-2))";
char Deriv_y_string[]                 = "kronDelta(1,var(-4,-1))*kronDelta(2,var(-4,-2))";
char Deriv_z_string[]                 = "kronDelta(1,var(-4,-1))*kronDelta(3,var(-4,-2))";
char Deriv_kronDelta_string[]         = "0";
char Deriv_diracDelta_string[]        = "var(-2,-1)*(diracDelta(var(-1,-1),var(-1,-2))/(var(-1,-1)-var(-1,-2))) - var(-2,-2)*(diracDelta(var(-1,-1),var(-1,-2))/(var(-1,-1)-var(-1,-2)))";
char Deriv_finDiracDel_string[]       = "var(-2,-1)*(finDiracDel(var(-1,-1),var(-1,-2))/(var(-1,-1)-var(-1,-2))) - var(-2,-2)*(finDiracDel(var(-1,-1),var(-1,-2))/(var(-1,-1)-var(-1,-2)))";
char Deriv_split_string[]             = "split(var(-1,-1),var(-1,-2),var(-2,-3),var(-2,-4))";
char Deriv_esplit_string[]            = "esplit(var(-1,-1),var(-1,-2),var(-2,-3),var(-2,-4))";
char Deriv_isint_string[]             = "isint(var(-1,-1),var(-2,-2),var(-2,-3))";
char Deriv_isreal_string[]            = "isreal(var(-1,-1),var(-2,-2),var(-2,-3))";
char Deriv_isanum_string[]            = "isanum(var(-1,-1),var(-2,-2),var(-2,-3))";
char Deriv_isfinint_string[]          = "isfinint(var(-1,-1),var(-2,-2),var(-2,-3))";
char Deriv_isfinreal_string[]         = "isfinreal(var(-1,-1),var(-2,-2),var(-2,-3))";
char Deriv_isfinanum_string[]         = "isfinanum(var(-1,-1),var(-2,-2),var(-2,-3))";
char Deriv_add_string[]               = "var(-2,-1)+var(-2,-2)";
char Deriv_mul_string[]               = "var(-2,-1)*var(-1,-2) + var(-1,-1)*var(-2,-2)";
char Deriv_div_string[]               = "var(-2,-1)/var(-1,-2) - var(-2,-2)*(var(-1,-1)/pow(var(-1,-2),2))";
char Deriv_idiv_string[]              = "0";
char Deriv_fmod_string[]              = "0";
char Deriv_eval_string[]              = "eval(var(-1,-1),var(-1,-2),var(-2,-3))";
char Deriv_sum_string[]               = "sum(var(-1,-1),var(-1,-2),var(-1,-3),var(-1,-4),var(-2,-5))";
char Deriv_prod_string[]              = "sum(var(-3,-1),var(-1,-2),var(-1,-3),var(-1,-4),prod(var(-1,-1),var(-1,-2),var(-1,-3),var(-1,-4),esplit(var(0,var(-3,-1)),var(0,var(-1,-1)),var(-2,-5),var(-1,-5))))";
char Deriv_integ_string[]             = "integ(var(-1,-1),var(-1,-2),var(-1,-3),var(-1,-4),var(-2,-5)) + (var(-2,-3)*eval(var(-1,-1),var(-1,-3),var(-1,-5))) - (var(-2,-2)*eval(var(-1,-1),var(-1,-2),var(-1,-5)))";
char Deriv_max_string[]               = "split(var(-1,-1),var(-1,-2),var(-2,-2),var(-2,-1))";
char Deriv_min_string[]               = "split(var(-1,-1),var(-1,-2),var(-2,-1),var(-2,-2))";
char Deriv_abs_string[]               = "var(-2,-1)*(diracDelta(var(-1,-1),0)+sgn(var(-1,-1)))";
char Deriv_sgn_string[]               = "var(-2,-1)*diracDelta(var(-1,-1),0)";
char Deriv_rint_string[]              = "var(-2,-1)*diracDelta(var(-1,-1),rint(var(-1,-1)))";
char Deriv_ceil_string[]              = "var(-2,-1)*diracDelta(var(-1,-1),ceil(var(-1,-1)))";
char Deriv_floor_string[]             = "var(-2,-1)*diracDelta(var(-1,-1),floor(var(-1,-1)))";
char Deriv_inv_string[]               = "var(-2,-1)*-pow(var(-1,-1),-2)";
char Deriv_neg_string[]               = "var(-2,-1)*-1";
char Deriv_sin_string[]               = "var(-2,-1)*cos(var(-1,-1))";
char Deriv_cos_string[]               = "var(-2,-1)*-sin(var(-1,-1))";
char Deriv_tan_string[]               = "var(-2,-1)*pow(sec(var(-1,-1)),2)";
char Deriv_cosec_string[]             = "var(-2,-1)*-cosec(var(-1,-1))*cot(var(-1,-1))";
char Deriv_sec_string[]               = "var(-2,-1)*sec(var(-1,-1))*tan(var(-1,-1))";
char Deriv_cot_string[]               = "var(-2,-1)*-pow(cosec(var(-1,-1)),2)";
char Deriv_asin_string[]              = "var(-2,-1)*inv(sqrt(1-pow(var(-1,-1),2)))";
char Deriv_acos_string[]              = "var(-2,-1)*-inv(sqrt(1-pow(var(-1,-1),2)))";
char Deriv_atan_string[]              = "var(-2,-1)*inv(1+pow(var(-1,-1),2))";
char Deriv_acosec_string[]            = "var(-2,-1)*-inv(var(-1,-1)*sqrt(pow(var(-1,-1),2)-1))";
char Deriv_asec_string[]              = "var(-2,-1)*inv(var(-1,-1)*sqrt(pow(var(-1,-1),2)-1))";
char Deriv_acot_string[]              = "var(-2,-1)*-inv(1+pow(var(-1,-1),2))";
char Deriv_sinh_string[]              = "var(-2,-1)*cosh(var(-1,-1))";
char Deriv_cosh_string[]              = "var(-2,-1)*sinh(var(-1,-1))";
char Deriv_tanh_string[]              = "var(-2,-1)*pow(sech(var(-1,-1)),2)";
char Deriv_cosech_string[]            = "var(-2,-1)*-cosech(var(-1,-1))*coth(var(-1,-1))";
char Deriv_sech_string[]              = "var(-2,-1)*-sech(var(-1,-1))*tanh(var(-1,-1))";
char Deriv_coth_string[]              = "var(-2,-1)*-pow(cosech(var(-1,-1)),2)";
char Deriv_asinh_string[]             = "var(-2,-1)*inv(sqrt(pow(var(-1,-1),2)+1))";
char Deriv_acosh_string[]             = "var(-2,-1)*inv(sqrt(pow(var(-1,-1),2)-1))";
char Deriv_atanh_string[]             = "var(-2,-1)*inv(1-pow(var(-1,-1),2))";
char Deriv_acosech_string[]           = "var(-2,-1)*-inv(var(-1,-1)*sqrt(pow(var(-1,-1),2)+1))";
char Deriv_asech_string[]             = "var(-2,-1)*-inv(var(-1,-1)*sqrt(1-pow(var(-1,-1),2)))";
char Deriv_acoth_string[]             = "var(-2,-1)*inv(1-pow(var(-1,-1),2))";
char Deriv_sinc_string[]              = "var(-2,-1)*((cos(3.14159*var(-1,-1))-sinc(var(-1,-1)))/var(-1,-1))";
char Deriv_gd_string[]                = "var(-2,-1)*sech(var(-1,-1))";
char Deriv_agd_string[]               = "var(-2,-1)*sec(var(-1,-1))";
char Deriv_exp_string[]               = "var(-2,-1)*exp(var(-1,-1))";
char Deriv_tenup_string[]             = "var(-2,-1)*log(10)*tenup(var(-1,-1))";
char Deriv_log_string[]               = "var(-2,-1)*inv(var(-1,-1))";
char Deriv_logX_string[]              = "var(-2,-1)*10*inv(var(-1,-1))";
char Deriv_pow_string[]               = "var(-2,-1)*var(-1,-2)*pow(var(-1,-1),var(-1,-2)-1) + var(-2,-2)*log(var(-1,-1))*pow(var(-1,-1),var(-1,-2))";
char Deriv_sqrt_string[]              = "var(-2,-1)*pow(var(-1,-1),-1/2)/2";
char Deriv_gamma_string[]             = "var(-2,-1)*gamma(var(-1,-1))*psi(var(-1,-1))";
char Deriv_lngamma_string[]           = "var(-2,-1)*psi(var(-1,-1))";
char Deriv_normDistr_string[]         = "var(-2,-1)*normDistr(var(-1,-1))*-var(-1,-1)";
char Deriv_polyDistr_string[]         = "var(-2,-2)*polyDistr(var(-1,-1),var(-1,-2))*-pow(var(-1,-2),var(-1,-1)-1)*pow(sgn(var(-1,-2)),var(-1,-1))*pow(sqrt(gamma(3/var(-1,-1))/gamma(1/var(-1,-1))),var(-1,-1))";
char Deriv_erf_string[]               = "var(-2,-1)*2*sqrt(2)*normDistr(sqrt(2)*var(-1,-1))";
char Deriv_erfc_string[]              = "var(-2,-1)*-2*sqrt(2)*normDistr(sqrt(2)*var(-1,-1))";
char Deriv_perm_string[]              = "0";
char Deriv_comb_string[]              = "0";
char Deriv_fact_string[]              = "0";
char Deriv_randn_string[]             = "0";
char Deriv_irand_string[]             = "0";
char Deriv_grand_string[]             = "0";
char Deriv_prand_string[]             = "0";
char Deriv_gami_string[]              = "var(-2,-1)*integ(var(-3,-1),0,var(-1,-2),0,log(var(0,var(-3,-1)))*pow(var(0,var(-3,-1)),var(-1,-1)-1)*exp(-var(0,var(-3,-1)))) + var(-2,-2)*pow(var(-1,-2),var(-1,-1)-1)*exp(-var(-1,-2))";
char Deriv_gamic_string[]             = "var(-2,-1)*(gamma(var(-1,-2))*psi(var(-1,-2)) - integ(var(-3,-1),0,var(-1,-2),0,log(var(0,var(-3,-1)))*pow(var(0,var(-3,-1)),var(-1,-1)-1)*exp(-var(0,var(-3,-1))))) - var(-2,-2)*pow(var(-1,-2),var(-1,-1)-1)*exp(-var(-1,-2))";
char Deriv_psi_string[]               = "var(-2,-1)*psi_n(1,var(-1,-1))";
char Deriv_psi_n_string[]             = "var(-2,-2)*psi_n(var(-1,-1)+1,var(-1,-2))";
char Deriv_dawson_string[]            = "var(-2,-1)*(1-(2*var(-1,-1)*dawson(var(-1,-1))))";

char Name_one[]         = "one";
char Name_zero[]        = "zero";
char Name_inf[]         = "inf";
char Name_ninf[]        = "ninf";
char Name_FiniteIndet[] = "FiniteIndet";
char Name_infintIndet[] = "infintIndet";
char Name_isError[]     = "isError";
char Name_var[]         = "var";
char Name_x[]           = "x";
char Name_y[]           = "y";
char Name_z[]           = "z";
char Name_randn[]       = "randn";
char Name_irand[]       = "irand";
char Name_grand[]       = "grand";
char Name_prand[]       = "prand";
char Name_split[]       = "split";
char Name_esplit[]      = "esplit";
char Name_isint[]       = "isint";
char Name_isreal[]      = "isreal";
char Name_isanum[]      = "isanum";
char Name_isfinint[]    = "isfinint";
char Name_isfinreal[]   = "isfinreal";
char Name_isfinanum[]   = "isfinanum";
char Name_max[]         = "max";
char Name_min[]         = "min";
char Name_inv[]         = "inv";
char Name_neg[]         = "neg";
char Name_add[]         = "add";
char Name_mul[]         = "mul";
char Name_div[]         = "div";
char Name_idiv[]        = "idiv";
char Name_fmod[]        = "fmod";
char Name_eval[]        = "eval";
char Name_sum[]         = "sum";
char Name_prod[]        = "prod";
char Name_integ[]       = "integ";
char Name_kronDelta[]   = "kronDelta";
char Name_diracDelta[]  = "diracDelta";
char Name_finDiracDel[] = "finDiracDel";
char Name_abs[]         = "abs";
char Name_sgn[]         = "sgn";
char Name_rint[]        = "rint";
char Name_ceil[]        = "ceil";
char Name_floor[]       = "floor";
char Name_perm[]        = "perm";
char Name_comb[]        = "comb";
char Name_fact[]        = "fact";
char Name_pow[]         = "pow";
char Name_sqrt[]        = "sqrt";
char Name_sin[]         = "sin";
char Name_cos[]         = "cos";
char Name_tan[]         = "tan";
char Name_cosec[]       = "cosec";
char Name_sec[]         = "sec";
char Name_cot[]         = "cot";
char Name_asin[]        = "asin";
char Name_acos[]        = "acos";
char Name_atan[]        = "atan";
char Name_acosec[]      = "acosec";
char Name_asec[]        = "asec";
char Name_acot[]        = "acot";
char Name_sinh[]        = "sinh";
char Name_cosh[]        = "cosh";
char Name_tanh[]        = "tanh";
char Name_cosech[]      = "cosech";
char Name_sech[]        = "sech";
char Name_coth[]        = "coth";
char Name_asinh[]       = "asinh";
char Name_acosh[]       = "acosh";
char Name_atanh[]       = "atanh";
char Name_acosech[]     = "acosech";
char Name_asech[]       = "asech";
char Name_acoth[]       = "acoth";
char Name_sinc[]        = "sinc";
char Name_gd[]          = "gd";
char Name_agd[]         = "agd";
char Name_exp[]         = "exp";
char Name_tenup[]       = "tenup";
char Name_log[]         = "log";
char Name_logX[]        = "logX";
char Name_gamma[]       = "gamma";
char Name_lngamma[]     = "lngamma";
char Name_gami[]        = "gami";
char Name_gamic[]       = "gamic";
char Name_psi[]         = "psi";
char Name_psi_n[]       = "psi_n";
char Name_erf[]         = "erf";
char Name_erfc[]        = "erfc";
char Name_dawson[]      = "dawson";
char Name_normDistr[]   = "normDistr";
char Name_polyDistr[]   = "polyDistr";








/*********************************************************************

Local Variables
===============

Local variables are dealt with using this set of functions.  As described
previously, local variables are things of the form var(0,i) where i>0.
There are MAX_NUM_LOC_VARS available.  Each is associated with two
variables, namely:

(*loc_vars)[i-1] = the contents of this local variables (uncalced default).
(*loc_var_used_flags)[i-1] = 0 normally (variable unused/free)
                             1 if variable has been set or touched.

All functions return non-zero on failure:

alloc_genarg_locals: allocate local variables in given context.  Returns
                     0 on success, -1 on failure.
dealloc_genarg_locals: deallocate the above.

find_first_unused_loc_var: find first unused local variable (ie. min i
                           such that loc_var_used_flags[i-1] = 0.
set_loc_var: set local variable.  Will also set loc_var_used_flags[i-1].
restore_loc_var: just like set_loc_var, but does not change the value
                 of loc_var_used_flags[i-1].
get_loc_var: get local variable (ignores flag).
touch_loc_var: sets flag for variable without changing value.

*********************************************************************/


int alloc_genarg_locals(GenArg *context);
void dealloc_genarg_locals(GenArg *context);
long find_first_unused_loc_var(GenArg context);
int set_loc_var(long which, GenVar what, GenArg context);
int restore_loc_var(long which, GenVar what, GenArg context);
GenVar get_loc_var(long which, GenArg context);
int touch_loc_var(long i, GenArg context);





/*********************************************************************

Miscellaneous
=============

The following functions provide random numbers as required.  Note that
I did not write these functions - attribution can be found later in the
code.

*********************************************************************/

void   RandomInitialise(int,int);
double RandomUniform(void);
double RandomGaussian(double mean, double stddev);
double RandomPolynomial(double mean, double stddev, int order);
int    RandomInt(int lower, int upper);
double RandomDouble(double lower, double upper);




/*********************************************************************

Function Definitions and LUT
============================

For every function *** there are two functions defined, namely an
evaluation function Calc_*** and a simplification function
Simplify_***.  Pointers to these are contained in the FunctionTable
LUT, in terms of an array of FunctionNodes.

*********************************************************************/

#define ONE_FN_NUMBER           0
#define ZERO_FN_NUMBER          1
#define PINF_FN_NUMBER          2
#define NINF_FN_NUMBER          3
#define FIND_FN_NUMBER          4
#define IIND_FN_NUMBER          5
#define ERR_FN_NUMBER           6

#define VAR_FN_NUMBER           7
#define X_FN_NUMBER             8
#define Y_FN_NUMBER             9
#define Z_FN_NUMBER             10

#define NUM_FUNCTIONS           91



int Calc_one(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_zero(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_inf(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);    
int Calc_ninf(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);    
int Calc_FiniteIndet(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly); 
int Calc_InfintIndet(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly); 
int Calc_isError(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);     
int Calc_var(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);
int Calc_x(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);
int Calc_y(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);
int Calc_z(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);
int Calc_kronDelta(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);   
int Calc_diracDelta(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);  
int Calc_finDiracDel(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);  
int Calc_split(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);       
int Calc_esplit(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);      
int Calc_isint(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);       
int Calc_isreal(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);      
int Calc_isanum(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);      
int Calc_isfinint(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);       
int Calc_isfinreal(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);      
int Calc_isfinanum(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);      
int Calc_add(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_mul(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_div(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_idiv(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_fmod(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_eval(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_sum(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_prod(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_integ(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);       
int Calc_max(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_min(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_abs(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_sgn(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_rint(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_ceil(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_floor(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);       
int Calc_inv(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_neg(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_sin(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_cos(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_tan(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_cosec(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);       
int Calc_sec(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_cot(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_asin(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_acos(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_atan(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_acosec(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);      
int Calc_asec(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_acot(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_sinh(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_cosh(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_tanh(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_cosech(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);      
int Calc_sech(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_coth(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_asinh(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);       
int Calc_acosh(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);       
int Calc_atanh(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);       
int Calc_acosech(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);     
int Calc_asech(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);       
int Calc_acoth(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);       
int Calc_sinc(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);          
int Calc_gd(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);          
int Calc_agd(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_exp(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_tenup(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);       
int Calc_log(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);          
int Calc_logX(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_pow(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_sqrt(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_gamma(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);       
int Calc_lngamma(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);      
int Calc_normDistr(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);   
int Calc_polyDistr(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);   
int Calc_erf(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_erfc(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_perm(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_comb(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);         
int Calc_fact(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);        
int Calc_randn(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);
int Calc_irand(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);       
int Calc_grand(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);
int Calc_prand(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);
int Calc_gami(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);
int Calc_gamic(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);
int Calc_psi(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);
int Calc_psi_n(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);
int Calc_dawson(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly);


FunctionNode FunctionTable[NUM_FUNCTIONS] = {

{ Name_one         , 0 , &Calc_one         , Deriv_one_string         , one_wderiv         , 0 , ',' } ,
{ Name_zero        , 0 , &Calc_zero        , Deriv_zero_string        , zero_wderiv        , 0 , ',' } ,
{ Name_inf         , 0 , &Calc_inf         , Deriv_inf_string         , inf_wderiv         , 0 , ',' } ,
{ Name_ninf        , 0 , &Calc_ninf        , Deriv_ninf_string        , ninf_wderiv        , 0 , ',' } ,
{ Name_FiniteIndet , 0 , &Calc_FiniteIndet , Deriv_FiniteIndet_string , FiniteIndet_wderiv , 0 , ',' } ,
{ Name_infintIndet , 0 , &Calc_InfintIndet , Deriv_InfintIndet_string , infintIndet_wderiv , 0 , ',' } ,
{ Name_isError     , 0 , &Calc_isError     , Deriv_isError_string     , isError_wderiv     , 0 , ',' } ,

{ Name_var         , 2 , &Calc_var         , Deriv_var_string         , var_wderiv         , 0 , ',' } ,
{ Name_x           , 0 , &Calc_x           , Deriv_x_string           , x_wderiv           , 0 , ',' } ,
{ Name_y           , 0 , &Calc_y           , Deriv_y_string           , y_wderiv           , 0 , ',' } ,
{ Name_z           , 0 , &Calc_z           , Deriv_z_string           , z_wderiv           , 0 , ',' } ,

{ Name_randn       , 2 , &Calc_randn       , Deriv_randn_string       , randn_wderiv       , 0 , ',' } ,
{ Name_irand       , 2 , &Calc_irand       , Deriv_irand_string       , irand_wderiv       , 0 , ',' } ,
{ Name_grand       , 2 , &Calc_grand       , Deriv_grand_string       , grand_wderiv       , 0 , ',' } ,
{ Name_prand       , 3 , &Calc_prand       , Deriv_prand_string       , prand_wderiv       , 0 , ',' } ,

{ Name_split       , 4 , &Calc_split       , Deriv_split_string       , split_wderiv       , 0 , ',' } ,
{ Name_esplit      , 4 , &Calc_esplit      , Deriv_esplit_string      , esplit_wderiv      , 0 , ',' } ,

{ Name_isint       , 3 , &Calc_isint       , Deriv_isint_string       , isint_wderiv       , 0 , ',' } ,
{ Name_isreal      , 3 , &Calc_isreal      , Deriv_isreal_string      , isreal_wderiv      , 0 , ',' } ,
{ Name_isanum      , 3 , &Calc_isanum      , Deriv_isanum_string      , isanum_wderiv      , 0 , ',' } ,
{ Name_isfinint    , 3 , &Calc_isfinint    , Deriv_isfinint_string    , isfinint_wderiv    , 0 , ',' } ,
{ Name_isfinreal   , 3 , &Calc_isfinreal   , Deriv_isfinreal_string   , isfinreal_wderiv   , 0 , ',' } ,
{ Name_isfinanum   , 3 , &Calc_isfinanum   , Deriv_isfinanum_string   , isfinanum_wderiv   , 0 , ',' } ,

{ Name_max         , 2 , &Calc_max         , Deriv_max_string         , max_wderiv         , 0 , ',' } ,
{ Name_min         , 2 , &Calc_min         , Deriv_min_string         , min_wderiv         , 0 , ',' } ,

{ Name_inv         , 1 , &Calc_inv         , Deriv_inv_string         , inv_wderiv         , 0 , ',' } ,
{ Name_neg         , 1 , &Calc_neg         , Deriv_neg_string         , neg_wderiv         , 0 , ',' } ,

{ Name_add         , 2 , &Calc_add         , Deriv_add_string         , add_wderiv         , 0 , '+' } ,
{ Name_mul         , 2 , &Calc_mul         , Deriv_mul_string         , mul_wderiv         , 0 , '*' } ,
{ Name_div         , 2 , &Calc_div         , Deriv_div_string         , div_wderiv         , 0 , '/' } ,
{ Name_idiv        , 2 , &Calc_idiv        , Deriv_idiv_string        , idiv_wderiv        , 0 , ',' } ,
{ Name_fmod        , 2 , &Calc_fmod        , Deriv_fmod_string        , fmod_wderiv        , 0 , ',' } ,

{ Name_eval        , 3 , &Calc_eval        , Deriv_eval_string        , eval_wderiv        , 0 , ',' } ,
{ Name_sum         , 5 , &Calc_sum         , Deriv_sum_string         , sum_wderiv         , 0 , ',' } ,
{ Name_prod        , 5 , &Calc_prod        , Deriv_prod_string        , prod_wderiv        , 1 , ',' } ,
{ Name_integ       , 5 , &Calc_integ       , Deriv_integ_string       , integ_wderiv       , 1 , ',' } ,

{ Name_kronDelta   , 2 , &Calc_kronDelta   , Deriv_kronDelta_string   , kronDelta_wderiv   , 0 , ',' } ,
{ Name_diracDelta  , 2 , &Calc_diracDelta  , Deriv_diracDelta_string  , diracDelta_wderiv  , 0 , ',' } ,
{ Name_finDiracDel , 2 , &Calc_finDiracDel , Deriv_finDiracDel_string , finDiracDel_wderiv , 0 , ',' } ,

{ Name_abs         , 1 , &Calc_abs         , Deriv_abs_string         , abs_wderiv         , 0 , ',' } ,
{ Name_sgn         , 1 , &Calc_sgn         , Deriv_sgn_string         , sgn_wderiv         , 0 , ',' } ,
{ Name_rint        , 1 , &Calc_rint        , Deriv_rint_string        , rint_wderiv        , 0 , ',' } ,
{ Name_ceil        , 1 , &Calc_ceil        , Deriv_ceil_string        , ceil_wderiv        , 0 , ',' } ,
{ Name_floor       , 1 , &Calc_floor       , Deriv_floor_string       , floor_wderiv       , 0 , ',' } ,

{ Name_perm        , 2 , &Calc_perm        , Deriv_perm_string        , perm_wderiv        , 0 , ',' } ,
{ Name_comb        , 2 , &Calc_comb        , Deriv_comb_string        , comb_wderiv        , 0 , ',' } ,
{ Name_fact        , 1 , &Calc_fact        , Deriv_fact_string        , fact_wderiv        , 0 , ',' } ,

{ Name_pow         , 2 , &Calc_pow         , Deriv_pow_string         , pow_wderiv         , 0 , '^' } ,
{ Name_sqrt        , 1 , &Calc_sqrt        , Deriv_sqrt_string        , sqrt_wderiv        , 0 , ',' } ,

{ Name_sin         , 1 , &Calc_sin         , Deriv_sin_string         , sin_wderiv         , 0 , ',' } ,
{ Name_cos         , 1 , &Calc_cos         , Deriv_cos_string         , cos_wderiv         , 0 , ',' } ,
{ Name_tan         , 1 , &Calc_tan         , Deriv_tan_string         , tan_wderiv         , 0 , ',' } ,
{ Name_cosec       , 1 , &Calc_cosec       , Deriv_cosec_string       , cosec_wderiv       , 0 , ',' } ,
{ Name_sec         , 1 , &Calc_sec         , Deriv_sec_string         , sec_wderiv         , 0 , ',' } ,
{ Name_cot         , 1 , &Calc_cot         , Deriv_cot_string         , cot_wderiv         , 0 , ',' } ,
{ Name_asin        , 1 , &Calc_asin        , Deriv_asin_string        , asin_wderiv        , 0 , ',' } ,
{ Name_acos        , 1 , &Calc_acos        , Deriv_acos_string        , acos_wderiv        , 0 , ',' } ,
{ Name_atan        , 1 , &Calc_atan        , Deriv_atan_string        , atan_wderiv        , 0 , ',' } ,
{ Name_acosec      , 1 , &Calc_acosec      , Deriv_acosec_string      , acosec_wderiv      , 0 , ',' } ,
{ Name_asec        , 1 , &Calc_asec        , Deriv_asec_string        , asec_wderiv        , 0 , ',' } ,
{ Name_acot        , 1 , &Calc_acot        , Deriv_acot_string        , acot_wderiv        , 0 , ',' } ,
{ Name_sinh        , 1 , &Calc_sinh        , Deriv_sinh_string        , sinh_wderiv        , 0 , ',' } ,
{ Name_cosh        , 1 , &Calc_cosh        , Deriv_cosh_string        , cosh_wderiv        , 0 , ',' } ,
{ Name_tanh        , 1 , &Calc_tanh        , Deriv_tanh_string        , tanh_wderiv        , 0 , ',' } ,
{ Name_cosech      , 1 , &Calc_cosech      , Deriv_cosech_string      , cosech_wderiv      , 0 , ',' } ,
{ Name_sech        , 1 , &Calc_sech        , Deriv_sech_string        , sech_wderiv        , 0 , ',' } ,
{ Name_coth        , 1 , &Calc_coth        , Deriv_coth_string        , coth_wderiv        , 0 , ',' } ,
{ Name_asinh       , 1 , &Calc_asinh       , Deriv_asinh_string       , asinh_wderiv       , 0 , ',' } ,
{ Name_acosh       , 1 , &Calc_acosh       , Deriv_acosh_string       , acosh_wderiv       , 0 , ',' } ,
{ Name_atanh       , 1 , &Calc_atanh       , Deriv_atanh_string       , atanh_wderiv       , 0 , ',' } ,
{ Name_acosech     , 1 , &Calc_acosech     , Deriv_acosech_string     , acosech_wderiv     , 0 , ',' } ,
{ Name_asech       , 1 , &Calc_asech       , Deriv_asech_string       , asech_wderiv       , 0 , ',' } ,
{ Name_acoth       , 1 , &Calc_acoth       , Deriv_acoth_string       , acoth_wderiv       , 0 , ',' } ,

{ Name_sinc        , 1 , &Calc_sinc        , Deriv_sinc_string        , sinc_wderiv        , 0 , ',' } , 
{ Name_gd          , 1 , &Calc_gd          , Deriv_gd_string          , gd_wderiv          , 0 , ',' } ,
{ Name_agd         , 1 , &Calc_agd         , Deriv_agd_string         , agd_wderiv         , 0 , ',' } ,

{ Name_exp         , 1 , &Calc_exp         , Deriv_exp_string         , exp_wderiv         , 0 , ',' } ,
{ Name_tenup       , 1 , &Calc_tenup       , Deriv_tenup_string       , tenup_wderiv       , 0 , ',' } ,
{ Name_log         , 1 , &Calc_log         , Deriv_log_string         , log_wderiv         , 0 , ',' } ,
{ Name_logX        , 1 , &Calc_logX        , Deriv_logX_string        , logX_wderiv        , 0 , ',' } ,

{ Name_gamma       , 1 , &Calc_gamma       , Deriv_gamma_string       , gamma_wderiv       , 0 , ',' } ,
{ Name_lngamma     , 1 , &Calc_lngamma     , Deriv_lngamma_string     , lngamma_wderiv     , 0 , ',' } ,
{ Name_gami        , 2 , &Calc_gami        , Deriv_gami_string        , gami_wderiv        , 1 , ',' } ,
{ Name_gamic       , 2 , &Calc_gamic       , Deriv_gamic_string       , gamic_wderiv       , 1 , ',' } ,
{ Name_psi         , 1 , &Calc_psi         , Deriv_psi_string         , psi_wderiv         , 0 , ',' } ,
{ Name_psi_n       , 2 , &Calc_psi_n       , Deriv_psi_n_string       , psi_n_wderiv       , 0 , ',' } ,

{ Name_erf         , 1 , &Calc_erf         , Deriv_erf_string         , erf_wderiv         , 0 , ',' } ,
{ Name_erfc        , 1 , &Calc_erfc        , Deriv_erfc_string        , erfc_wderiv        , 0 , ',' } ,
{ Name_dawson      , 1 , &Calc_dawson      , Deriv_dawson_string      , dawson_wderiv      , 0 , ',' } ,
{ Name_normDistr   , 1 , &Calc_normDistr   , Deriv_normDistr_string   , normDistr_wderiv   , 0 , ',' } ,
{ Name_polyDistr   , 2 , &Calc_polyDistr   , Deriv_polyDistr_string   , polyDistr_wderiv   , 0 , ',' }

};


















































/*********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
*****************************      ***********************************
************************                ******************************
*********************                      ***************************
********************    Code starts here    **************************
********************    ================    **************************
*********************                      ***************************
************************                ******************************
*****************************      ***********************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
*********************************************************************/






/*********************************************************************

   Functions: parsing, deparsing, copying and deleting
   ===================================================

*********************************************************************/

MathsNode *parseMaths(char *express)
{
    MathsNode *result;

    if ( express == NULL )
    {
        return NULL;
    }

    if ( ( result = fast_parseMaths(express) ) != NULL )
    {
        if ( simplifyEqn(result) )
        {
            delEqn(result);

            return NULL;
        }
    }

    return result;
}

MathsNode *fast_parseMaths(char *express)
{
    MathsNode *result;
    char *form_express[1];
    long len;

    if ( express == NULL )
    {
        return NULL;
    }

    if ( ( len = convert_equation(express,form_express) ) <= 0 )
    {
        return NULL;
    }

    result = locparseMaths(*form_express,len);

    /*
       Errors here will only mean result is NULL, but won't affect
       form_express.
    */

    free(*form_express);

    return result;
}

MathsNode *locparseMaths(char *express, long len)
{
    long i,j,k,l;
    unsigned long fn_length;
    unsigned long mul_start;
    unsigned long mul_length;
    unsigned long arg_start;
    unsigned long arg_length;
    MathsNode *result;
    GenVar confuse;

    /*
       Assuming that the function has form FnName[_1_](_2_)..., where len is
       the length of FnName[_1_](_2_) (no null termination here) find the #
       of character before [ (i) and the number of chars before ( (j).
    */

    i = len_to_open(express,len);
    j = len_to_br_open(express,len);

    if ( i <= 0   ) { return NULL; }
    if ( j <= 0   ) { return NULL; }

    if ( j   < i+2 ) { return NULL; }
    if ( len < j+2 ) { return NULL; }

    /*
       fn_length        = # chars in FnName (must be >= 1)
       FnName+mul_start = "_1_](_2_)..."
       mul_length       = # chars in _1_ (must be >= 1)
       FnName+arg_start = "_2_)..."
       arg_length       = # chars in _2_ (must be >= 0)
    */

    fn_length  = i;
    mul_start  = i+1;
    mul_length = j-(i+2);
    arg_start  = j+1;
    arg_length = len-(j+2);

    if ( fn_length  < 1 ) { return NULL; }
    if ( mul_length < 1 ) { return NULL; }

    if ( ( result = (MathsNode *) malloc(sizeof(MathsNode)) ) == NULL )
    {
         return NULL;
    }

    /*
       First step - find what function this is.  -1 is the fall through
       case if the function is not recognised.
    */

    result->CombinerFnNumber = -1;

    for ( i = 0 ; i < NUM_FUNCTIONS ; i++ )
    {
        if ( ( strncmp(express,((FunctionTable[i]).FnName),fn_length) == 0 ) &&
             ( fn_length == strlen((FunctionTable[i]).FnName) ) )
        {
            result->CombinerFnNumber = i;

            break;
        }
    }

    if ( result->CombinerFnNumber == -1 )
    {
        free(result);

        return NULL;
    }

    result->multiplier = parseNumber(express+mul_start,mul_length);

    /*
       If the number is a constant, the multiplier must be subsumed into
       the number itself.  For example, if the multiplier is inf and the
       function is one, then the function must be changed to inf.
    */

    switch ( result->CombinerFnNumber )
    {
        case ONE_FN_NUMBER:
        case ZERO_FN_NUMBER:
        case PINF_FN_NUMBER:
        case NINF_FN_NUMBER:
        case FIND_FN_NUMBER:
        case IIND_FN_NUMBER:
        {
            confuse = ResultOne();

            make_number(result,&confuse);

            break;
        }

        default:
        {
            break;
        }
    }

    /*
       No work down the tree to allocate the subfunctions.
    */

    result->SubFunctions = NULL;

    j = (FunctionTable[result->CombinerFnNumber]).NumLocInput;

    if ( j > 0 )
    {
        if ( ( result->SubFunctions = (MathsNode **) malloc(j*sizeof(MathsNode *)) ) == NULL )
        {
            free(result);

            return NULL;
        }

        k = arg_start;

        for ( i = 0 ; i < j ; i++ )
        {
            /*
               Find the length of this argument, which starts are express+k
               and has max length len-k.  This is done by finding the length
               to the next , if it isn't the final argument, or the next )
               otherwise.
            */

            if ( i < j-1 )
            {
                l = len_to_comma(express+k,len-k);
            }

            else
            {
                l = len_to_br_close(express+k,len-k);

                /*
                   Sanity check - we should not be able to find any commas
                   before this bracket.  If we can, then there are too many
                   arguments in this function, so get the hell out of here.
                */

                if ( len_to_comma(express+k,len-k) > 0 )
                {
                    if ( i > 0 )
                    {
                        for ( j = 0 ; j < i ; j++ )
                        {
                            free((result->SubFunctions)[j]);
                        }
                    }

                    free(result->SubFunctions);
                    free(result);

                    return NULL;
                }
            }

            /*
               If something bad has happened, leave.
            */

            if ( l <= 0 )
            {
                if ( i > 0 )
                {
                    for ( j = 0 ; j < i ; j++ )
                    {
                        free((result->SubFunctions)[j]);
                    }
                }

                free(result->SubFunctions);
                free(result);

                return NULL;
            }

            /*
               OK, everything seems ok.  This subfunction should start
               at express+k and have length l.  Try to evaluate it.  If
               it fails, take the usual evasive action.
            */

            if ( ( (result->SubFunctions)[i] = locparseMaths(express+k,l) ) == NULL )
            {
                if ( i > 0 )
                {
                    for ( j = 0 ; j < i ; j++ )
                    {
                        free((result->SubFunctions)[j]);
                    }
                }

                free(result->SubFunctions);
                free(result);

                return NULL;
            }

            /*
               Update counters so that express+k refers to the start of the
               next argument.
            */

            k += l;
            k++;
        }

        /*
           If everything has gone ok then k should be at the end of the
           function.  If it isn't, report an error.
        */

        if ( k != len )
        {
            for ( i = 0 ; i < j ; i++ )
            {
                free((result->SubFunctions)[i]);
            }

            free(result->SubFunctions);
            free(result);

            return NULL;
        }
    }

    return result;
}

MathsNode *getZero(void)
{
    return copyEquation(xglobal_zero);
}

char *deparseMaths(MathsNode *equation)
{
    long i,j;
    long len;
    char *result;
    char *charmul;
    char **subresult;
    int neg_switch = 0;

    if ( equation == NULL )
    {
        return NULL;
    }

    /*
       What get's printed depends on what the function is.  For human
       readability, ONE_FN_NUMBER is printed as 1, and ZERO_FN_NUMBER
       as 0.
    */

    switch ( equation->CombinerFnNumber )
    {
        case ONE_FN_NUMBER:
        {
            if ( ( result = deparseNumber(equation->multiplier) ) == NULL )
            {
                return NULL;
            }

            break;
        }

        case ZERO_FN_NUMBER:
        {
            if ( ( result = (char *) malloc(2*sizeof(char)) ) == NULL )
            {
                return NULL;
            }

            result[0] = '0';
            result[1] = '\0';

            break;
        }

        default:
        {
            /*
               OK, the function is non-trivial (kinda).  What happens here
               really depends on what separator is being used.  If it is a
               comma then the function will be printed in standard form (but
               no [], as that's not human readable).  Otherwise, things are
               different.
            */

            switch ( (FunctionTable[equation->CombinerFnNumber]).altMarker )
            {
                case ',':
                {
                    /*
                       Standard function deparsing.  result is:
                       _0_*FnName(_1_,_2_,...)
                       where _0_ is the multiplier (presumed to be a number)
                       and _1_,_2_,... are the arguments.

                       len = the predicted length of the string.
                    */

                    len = strlen((FunctionTable[equation->CombinerFnNumber]).FnName);

                    j = (FunctionTable[equation->CombinerFnNumber]).NumLocInput;

                    if ( j > 0 )
                    {
                        /*
                           There are arguments - we need to deparse these
                           before we can do anything.
                        */

                        if ( ( subresult = (char **) malloc(j*sizeof(char *)) ) == NULL )
                        {
                            return NULL;
                        }

                        for ( i = 0 ; i < j ; i++ )
                        {
                            if ( ( subresult[i] = deparseMaths((equation->SubFunctions)[i]) ) == NULL )
                            {
                                if ( i > 0 )
                                {
                                    for ( j = 0 ; j < i ; j++ )
                                    {
                                        free(subresult[j]);
                                    }
                                }

                                free(subresult);

                                return NULL;
                            }

                            /*
                               The length is increased by the length of this
                               argument plus 1 for a ,.
                            */

                            len += strlen(subresult[i]);
                            len++;
                        }

                        /*
                           Take 1 from length (no , after the last arg).
                        */

                        len--;

                        /*
                           Brackets take 2 chars, so increase len.
                        */

                        len += 2;
                    }

                    else
                    {
                        subresult = NULL;
                    }

                    /*
                       For readability, we don't want to have something
                       like 1*sin(x), so we don't bother with a multiplier
                       if it is 1.
                    */

                    charmul = NULL;

                    if ( (equation->multiplier).DataType != DATA_IS_ONE )
                    {
                        if ( ( (equation->multiplier).DataType == DATA_IS_INTEGER ) &&
                             ( (equation->multiplier).IntegerValue == -1          )    )
                        {
                            /*
                               We need to replace FnName(...) with
                               -FnName(...), so increase len to suit.  Also
                               set flag.
                            */

                            neg_switch = 1;

                            len++;
                        }

                        else
                        {
                            /*
                               Having a non-trivial multiplier, we need to
                               convert it to a string.
                            */

                            if ( ( charmul = deparseNumber(equation->multiplier) ) == NULL )
                            {
                                if ( j > 0 )
                                {
                                    for ( i = 0 ; i < j ; i++ )
                                    {
                                        free(subresult[i]);
                                    }

                                    free(subresult);
                                }

                                return NULL;
                            }

                            /*
                               We need to replace FnName(...) with
                               (_0_*FnName(...)), so increase len to suit.
                            */

                            len += strlen(charmul);
                            len += 3;
                        }
                    }

                    if ( ( result = (char *) malloc((len+1)*sizeof(char)) ) == NULL )
                    {
                        if ( j > 0 )
                        {
                            for ( i = 0 ; i < j ; i++ )
                            {
                                free(subresult[i]);
                            }

                            free(subresult);
                        }

                        if ( charmul != NULL )
                        {
                            free(charmul);
                        }

                        return NULL;
                    }

                    if ( j > 0 )
                    {
                        /*
                           OK, we have all we need - now put it together.
                           First, the function, along with multipliers and
                           brackets.
                        */

                        if ( charmul == NULL )
                        {
                            if ( neg_switch )
                            {
                                sprintf(result,"-%s(",(FunctionTable[equation->CombinerFnNumber]).FnName);
                            }

                            else
                            {
                                sprintf(result,"%s(",(FunctionTable[equation->CombinerFnNumber]).FnName);
                            }
                        }

                        else
                        {
                            sprintf(result,"(%s*%s(",charmul,(FunctionTable[equation->CombinerFnNumber]).FnName);
                        }

                        /*
                           Then arguments.  Note that the first one does not
                           get preceeded by a comma.
                        */

                        for ( i = 0 ; i < j ; i++ )
                        {
                            if ( i )
                            {
                                sprintf(result,"%s%c%s",result,(FunctionTable[equation->CombinerFnNumber]).altMarker,subresult[i]);
                            }

                            else
                            {
                                sprintf(result,"%s%s",result,subresult[i]);
                            }
                        }

                        /*
                           Then finish it off.
                        */

                        sprintf(result,"%s)",result);

                        if ( charmul != NULL )
                        {
                            sprintf(result,"%s)",result);
                        }

                        /*
                            and free memory.
                        */

                        for ( i = 0 ; i < j ; i++ )
                        {
                            free(subresult[i]);
                        }

                        free(subresult);

                        if ( charmul != NULL )
                        {
                            free(charmul);
                        }
                    }

                    else
                    {
                        if ( charmul == NULL )
                        {
                            if ( neg_switch )
                            {
                                sprintf(result,"-%s",(FunctionTable[equation->CombinerFnNumber]).FnName);
                            }

                            else
                            {
                                sprintf(result,"%s",(FunctionTable[equation->CombinerFnNumber]).FnName);
                            }
                        }

                        else
                        {
                            sprintf(result,"(%s*%s)",charmul,(FunctionTable[equation->CombinerFnNumber]).FnName);

                            free(charmul);
                        }
                    }

                    break;
                }

                default:
                {
                    len = 0;

                    j = (FunctionTable[equation->CombinerFnNumber]).NumLocInput;

                    if ( j > 0 )
                    {
                        if ( ( subresult = (char **) malloc(j*sizeof(char *)) ) == NULL )
                        {
                            return NULL;
                        }

                        for ( i = 0 ; i < j ; i++ )
                        {
                            if ( ( subresult[i] = deparseMaths((equation->SubFunctions)[i]) ) == NULL )
                            {
                                if ( i > 0 )
                                {
                                    for ( j = 0 ; j < i ; j++ )
                                    {
                                        free(subresult[j]);
                                    }
                                }

                                free(subresult);

                                return NULL;
                            }

                            len += strlen(subresult[i]);
                            len++;
                        }

                        len--;

                        len += 2;
                    }

                    else
                    {
                        subresult = NULL;
                    }

                    charmul = NULL;

                    if ( (equation->multiplier).DataType != DATA_IS_ONE )
                    {
                        if ( ( (equation->multiplier).DataType == DATA_IS_INTEGER ) &&
                             ( (equation->multiplier).IntegerValue == -1          )    )
                        {
                            /*
                               We need to replace FnName(...) with
                               -FnName(...), so increase len to suit.  Also
                               set flag.
                            */

                            neg_switch = 1;

                            len++;
                        }

                        else
                        {
                            if ( ( charmul = deparseNumber(equation->multiplier) ) == NULL )
                            {
                                if ( j > 0 )
                                {
                                    for ( i = 0 ; i < j ; i++ )
                                    {
                                        free(subresult[i]);
                                    }

                                    free(subresult);
                                }

                                return NULL;
                            }

                           len += strlen(charmul);
                           len += 3;
                        }
                    }

                    if ( ( result = (char *) malloc((len+1)*sizeof(char)) ) == NULL )
                    {
                        if ( j > 0 )
                        {
                            for ( i = 0 ; i < j ; i++ )
                            {
                                free(subresult[i]);
                            }

                            free(subresult);
                        }

                        if ( charmul != NULL )
                        {
                            free(charmul);
                        }

                        return NULL;
                    }

                    if ( j > 0 )
                    {
                        if ( charmul == NULL )
                        {
                            if ( neg_switch )
                            {
                                sprintf(result,"-(");
                            }

                            else
                            {
                                sprintf(result,"(");
                            }
                        }

                        else
                        {
                            sprintf(result,"(%s*(",charmul);
                        }

                        for ( i = 0 ; i < j ; i++ )
                        {
                            if ( i )
                            {
                                sprintf(result,"%s%c%s",result,(FunctionTable[equation->CombinerFnNumber]).altMarker,subresult[i]);
                            }

                            else
                            {
                                sprintf(result,"%s%s",result,subresult[i]);
                            }
                        }

                        sprintf(result,"%s)",result);

                        for ( i = 0 ; i < j ; i++ )
                        {
                            free(subresult[i]);
                        }

                        free(subresult);

                        if ( charmul != NULL )
                        {
                            sprintf(result,"%s)",result);

                            free(charmul);
                        }
                    }

                    else
                    {
                        if ( j > 0 )
                        {
                            for ( i = 0 ; i < j ; i++ )
                            {
                                free(subresult[i]);
                            }

                            free(subresult);
                        }

                        if ( charmul != NULL )
                        {
                            free(charmul);
                        }

                        return NULL;
                    }

                    break;
                }
            }

            break;
        }
    }

    return result;
}

char *deparseMaths_cstruct(MathsNode *equation, char *name)
{
    int namelen;
    int len;
    int intlen;
    int intexp;
    char *tempa;
    char *result;
    int index = 0;

    if ( ( equation == NULL ) || ( name == NULL ) )
    {
        return NULL;
    }

    /*
       Stategy: the function locdeparseMaths_cstruct will create a tree
       of local variables, the top level being "MathsNode name_content_*_,
       where * is the value of index-1 (the index is used to create local
       variables in a valid namespace).  We then create a pointer to this
       with the correct name.

       NOTE: in this context, the variable index is a quasi-global counter.
             Because we don't know the structure of the equation except
             locally at any one time, we need a global counter to prevent
             variable name clashes.  So, whenever a variable is named, the
             value of index is incorporated into the name and index is
             incremented.  Consequently, we can work out the name of the
             variable allocated by locdeparseMaths_cstruct (at highest
             level) using index and our knowledge of the code.  And no, I
             can't be bothered making that explanation make strict sense.
    */

    namelen = strlen(name);

    if ( ( tempa = locdeparseMaths_cstruct(equation,name,namelen,&index) ) == NULL )
    {
        return NULL;
    }

    /*
       Calculate how many characters are needed to hold the string form
       of index, and store the result in intlen.
    */

    intlen = 1;
    intexp = 10;

    repeat_point:

    if ( index >= intexp )
    {
        intlen++;
        intexp *= 10;

        goto repeat_point;
    }

    /*
       Calculate the length of the result string (including termination and
       newline).  This depends on the sprintf's below.  Allow that \n may
       not be 1 character in dos (I really don't know if this matters -
       after all, '\n' is perfectly valid, and that *must* be one char).

       The magic number 30 is the amount of space taken by other "stuff".
    */

    len = 30 + strlen(tempa) + 2*namelen + intlen;

    /*
       Allocate and then create the string.
    */

    if ( ( result = (char *) malloc((len+1)*sizeof(char)) ) == NULL )
    {
        free(tempa);

        return NULL;
    }

    sprintf(result,"%s",tempa);
    sprintf(result,"%sMathsNode *",result);
    sprintf(result,"%s%s",result,name);
    sprintf(result,"%s = &",result);
    sprintf(result,"%s%s",result,name);
    sprintf(result,"%s_content_",result);
    sprintf(result,"%s%d",result,index-1);
    sprintf(result,"%s_;\n",result);

    free(tempa);

    return result;
}

char *locdeparseMaths_cstruct(MathsNode *equation, char *name, int namelen, int *index)
{
    char *result;
    char *temp;
    int num_subs;
    int i,j;
    int len;
    int indexlen = 1;
    int indexexp = 10;
    int fnnlen = 1;
    int fnnexp = 10;
    int sublen = 1;
    int subexp = 10;
    char **local_result;
    int *local_index;

    num_subs = (FunctionTable[equation->CombinerFnNumber]).NumLocInput;

    /*
       We need to know the length of all relevant integers, so we know how
       much space to allocate in the string.
    */

    repeat_point_a:

    if ( (*index) >= indexexp )
    {
        indexlen++;
        indexexp *= 10;

        goto repeat_point_a;
    }

    repeat_point_b:

    if ( (equation->CombinerFnNumber) >= fnnexp )
    {
        fnnlen++;
        fnnexp *= 10;

        goto repeat_point_b;
    }

    repeat_point_c:

    if ( num_subs >= subexp )
    {
        sublen++;
        subexp *= 10;

        goto repeat_point_c;
    }

    /*
       Get a string version of the multiplier.
    */

    if ( ( temp = deparseNumber_cstruct(equation->multiplier) ) == NULL )
    {
        return NULL;
    }

    if ( num_subs == 0 )
    {
        /*
           If there are no subfunctions then this is trivial.  The result
           is just "MathsNode *name_content_*_ = { * , * , NULL };", where
           the *'s are filled in appropriately and name is replaced with
           the name string.
        */

        len = ( 40 + namelen + indexlen + fnnlen + strlen(temp) );

        if ( ( result = (char *) malloc((len+1)*sizeof(char)) ) == NULL )
        {
            free(temp);

            return NULL;
        }

        sprintf(result,"MathsNode %s_content_%d_ = { %d , %s , NULL };\n",name,*index,equation->CombinerFnNumber,temp);
    }

    else
    {
        /*
           If there are subfunctions, these need to be statically allocated
           first (in the c code, above the current definition) so we have
           something we can point to (wtf can't c deal with pointers into
           the future?).  So lets get all relevant data, including the
           indexes.
        */

        if ( ( local_result = (char **) malloc(num_subs*sizeof(char *)) ) == NULL )
        {
            free(temp);

            return NULL;
        }

        if ( ( local_index = (int *) malloc(num_subs*sizeof(int)   ) ) == NULL )
        {
            free(local_result);
            free(temp);

            return NULL;
        }

        for ( i = 1 ; i <= num_subs ; i++ )
        {
            if ( ( local_result[i-1] = locdeparseMaths_cstruct((equation->SubFunctions)[i-1],name,namelen,index) ) == NULL )
            {
                if ( i > 1 )
                {
                    for ( j = 1 ; j < i ; j++ )
                    {
                        free(local_result[i-1]);
                    }
                }

                free(local_result);
                free(local_index);
                free(temp);

                return NULL;
            }

            local_index[i-1] = (*index)-1;
        }

        /*
           The index length will have changed, so we need to re-calculate
           the length of the index string.
        */

        indexlen = 1;
        indexexp = 10;

        repeat_point_ab:

        if ( (*index) >= indexexp )
        {
            indexlen++;
            indexexp *= 10;

            goto repeat_point_ab;
        }

        /*
           pre-calculate the length of the string.
        */

        len  = ( 28 + namelen + indexlen + sublen );
        len += ( 44 + namelen + indexlen + fnnlen + strlen(temp) + namelen + indexlen );

        for ( i = 1 ; i <= num_subs ; i++ )
        {
            len += strlen(local_result[i-1]);
            len += ( 15 + namelen + indexlen + sublen );
        }

        if ( ( result = (char *) malloc((len+1)*sizeof(char)) ) == NULL )
        {
            for ( j = 1 ; j < i ; j++ )
            {
                free(local_result[i-1]);
            }

            free(local_result);
            free(local_index);
            free(temp);

            return NULL;
        }

        for ( i = 1 ; i <= num_subs ; i++ )
        {
            if ( i == 1 )
            {
                sprintf(result,"%s",local_result[i-1]);
            }

            else
            {
                sprintf(result,"%s%s",result,local_result[i-1]);
            }
        }

        sprintf(result,"%sMathsNode *%s_refer_%d_[%d] = { ",result,name,*index,num_subs);

        for ( i = 1 ; i <= num_subs ; i++ )
        {
            if ( i < num_subs )
            {
                sprintf(result,"%s&%s_content_%d_ , ",result,name,local_index[i-1]);
            }

            else
            {
                sprintf(result,"%s&%s_content_%d_ };\n",result,name,local_index[i-1]);
            }
        }

        sprintf(result,"%sMathsNode %s_content_%d_ = { %d , %s , %s_refer_%d_ };\n",result,name,*index,equation->CombinerFnNumber,temp,name,*index);

        for ( i = 1 ; i <= num_subs ; i++ )
        {
            free(local_result[i-1]);
        }

        free(local_result);
        free(local_index);
    }

    free(temp);

    /*
       Make sure that there are no namespace clashes.
    */

    (*index)++;

    return result;
}

MathsNode *copyEquation(MathsNode *source)
{
    long i,j;
    MathsNode *result;

    if ( source == NULL )
    {
        return NULL;
    }

    if ( ( result = (MathsNode *) malloc(sizeof(MathsNode)) ) == NULL )
    {
        return NULL;
    }

    result->CombinerFnNumber = source->CombinerFnNumber;
    result->multiplier       = source->multiplier;

    result->SubFunctions = NULL;

    j = (FunctionTable[result->CombinerFnNumber]).NumLocInput;

    if ( j > 0 )
    {
        if ( ( result->SubFunctions = (MathsNode **) malloc(j*sizeof(MathsNode *)) ) == NULL )
        {
            free(result);

            return NULL;
        }

        for ( i = 0 ; i < j ; i++ )
        {
            if ( ( (result->SubFunctions)[i] = copyEquation((source->SubFunctions)[i]) ) == NULL )
            {
                if ( i > 0 )
                {
                    for ( j = 0 ; j < i ; j++ )
                    {
                        free((result->SubFunctions)[j]);
                    }
                }

                free(result->SubFunctions);
                free(result);

                return NULL;
            }
        }
    }

    return result;
}

int delEqn(MathsNode *equation)
{
    long i,j;

    if ( equation != NULL )
    {
        j = (FunctionTable[equation->CombinerFnNumber]).NumLocInput;

        if ( j > 0 )
        {
            if ( equation->SubFunctions != NULL )
            {
                for ( i = 0 ; i < j ; i++ )
                {
                    if ( (equation->SubFunctions)[i] != NULL )
                    {
                        delEqn((equation->SubFunctions)[i]);
                    }
                }

                free(equation->SubFunctions);
            }
        }

        free(equation);
    }

    return 0;
}




/*********************************************************************

   Functions: manipulation
   =======================

*********************************************************************/

int simplifyEqn(MathsNode *equation)
{
    GenVar result;
    int q;

    /*
       The NULL evaluation function does our work here - we just need
       to tell it to simplify on the fly and not to evaluate random
       numbers (just return DATA_IS_UNCALCED).
    */

    if ( ( q = Gen_evaluateEqnNull(&result,equation,1,1) ) )
    {
        return q;
    }

    return 0;
}

int subVar(MathsNode *equation, MathsNode *var_replace, long i, long j, int naive)
{
    int q;

    if ( ( equation == NULL ) || ( var_replace == NULL ) )
    {
        return EERR_BADARGS;
    }

    if ( ( i < 0 ) || ( j < 0 ) )
    {
        return EERR_BADARGS;
    }

    if ( ( q = locsubVar(equation,var_replace,i,j,naive) ) )
    {
        return q;
    }

    if ( ( q = simplifyEqn(equation) ) )
    {
        return q;
    }

    return 0;
}

int subVarRow(MathsNode *equation, MathsNode *var_replace, long i, int naive)
{
    int q;

    if ( ( equation == NULL ) || ( var_replace == NULL ) )
    {
        return EERR_BADARGS;
    }

    if ( i < 0 )
    {
        return EERR_BADARGS;
    }

    if ( ( q = locsubVarRow(equation,var_replace,i,naive) ) )
    {
        return q;
    }

    if ( ( q = simplifyEqn(equation) ) )
    {
        return q;
    }

    return 0;
}

int subVarCol(MathsNode *equation, MathsNode *var_replace, long j, int naive)
{
    int q;

    if ( ( equation == NULL ) || ( var_replace == NULL ) )
    {
        return EERR_BADARGS;
    }

    if ( j < 0 )
    {
        return EERR_BADARGS;
    }

    if ( ( q = locsubVarCol(equation,var_replace,j,naive) ) )
    {
        return q;
    }

    if ( ( q = simplifyEqn(equation) ) )
    {
        return q;
    }

    return 0;
}

int fast_subVar(MathsNode *equation, MathsNode *var_replace, long i, long j, int naive)
{
    if ( ( equation == NULL ) || ( var_replace == NULL ) )
    {
        return EERR_BADARGS;
    }

    if ( ( i < 0 ) || ( j < 0 ) )
    {
        return EERR_BADARGS;
    }

    return locsubVar(equation,var_replace,i,j,naive);
}

int fast_subVarRow(MathsNode *equation, MathsNode *var_replace, long i, int naive)
{
    if ( ( equation == NULL ) || ( var_replace == NULL ) )
    {
        return EERR_BADARGS;
    }

    if ( i < 0 )
    {
        return EERR_BADARGS;
    }

    return locsubVarRow(equation,var_replace,i,naive);
}

int fast_subVarCol(MathsNode *equation, MathsNode *var_replace, long j, int naive)
{
    if ( ( equation == NULL ) || ( var_replace == NULL ) )
    {
        return EERR_BADARGS;
    }

    if ( j < 0 )
    {
        return EERR_BADARGS;
    }

    return locsubVarCol(equation,var_replace,j,naive);
}

int locsubVar(MathsNode *equation, MathsNode *var_replace, long i, long j, int naive)
{
    long k,l;
    int q;
    GenVar tempi;
    GenVar tempj;
    MathsNode *resulta;
    MathsNode *resultb;
    MathsNode **tempz;

    k = (FunctionTable[equation->CombinerFnNumber]).NumLocInput;

    /*
       First substitute in any arguments.
    */

    if ( k > 0 )
    {
        for ( l = 0 ; l < k ; l++ )
        {
            if ( ( q = locsubVar((equation->SubFunctions)[l],var_replace,i,j,naive) ) )
            {
                return q;
            }
        }
    }

    if ( ( equation->CombinerFnNumber == VAR_FN_NUMBER ) ||
         ( equation->CombinerFnNumber == X_FN_NUMBER   ) ||
         ( equation->CombinerFnNumber == Y_FN_NUMBER   ) ||
         ( equation->CombinerFnNumber == Z_FN_NUMBER   )    )
    {
        /*
           OK, we've hit a variable of some sort.  Find out which one.
        */

        switch ( equation->CombinerFnNumber )
        {
            case X_FN_NUMBER:
            {
                tempi = ResultOne();
                tempj = ResultOne();

                equation->SubFunctions = NULL;

                break;
            }

            case Y_FN_NUMBER:
            {
                tempi = ResultOne();
                tempj = ResultInteger(2);

                equation->SubFunctions = NULL;

                break;
            }

            case Z_FN_NUMBER:
            {
                tempi = ResultOne();
                tempj = ResultInteger(3);

                equation->SubFunctions = NULL;

                break;
            }

            default:
            {
                if ( ( q = Gen_evaluateEqnNull(&tempi,(equation->SubFunctions)[0],1,0) ) )
                {
                    return q;
                }

                if ( ( q = Gen_evaluateEqnNull(&tempj,(equation->SubFunctions)[1],1,0) ) )
                {
                    return q;
                }

                break;
            }
        }

        if ( isFiniteInteger(tempi) && isFiniteInteger(tempj) )
        {
            /*
               Easiest case - it's a variable where both arguments are
               integers (ie. known) - in this case it is straightforward
               to replace it.
            */

            if ( ( getInteger(tempi) == i ) && ( getInteger(tempj) == j ) )
            {
                /*
                   Save the current variable for future use.
                */

                tempz = equation->SubFunctions;

                equation->CombinerFnNumber = var_replace->CombinerFnNumber;
                equation->multiplier       = mulGenVar(equation->multiplier,var_replace->multiplier);

                k = (FunctionTable[equation->CombinerFnNumber]).NumLocInput;

                if ( k > 0 )
                {
                    if ( ( equation->SubFunctions = (MathsNode **) malloc(k*sizeof(MathsNode *)) ) == NULL )
                    {
                        equation->SubFunctions = tempz;

                        return EERR_MALLOCFAIL;
                    }

                    for ( l = 0 ; l < k ; l++ )
                    {
                        if ( ( (equation->SubFunctions)[l] = copyEquation((var_replace->SubFunctions)[l]) ) == NULL )
                        {
                            if ( l > 0 )
                            {
                                for ( k = 0 ; k < l ; k++ )
                                {
                                    delEqn((equation->SubFunctions)[k]);
                                }
                            }

                            free(equation->SubFunctions);

                            equation->SubFunctions = tempz;

                            return EERR_BADCOPY;
                        }
                    }
                }

                else
                {
                    equation->SubFunctions = NULL;
                }

                /*
                   Delete old variable.
                */

                if ( tempz != NULL )
                {
                    delEqn(tempz[0]);
                    delEqn(tempz[1]);

                    free(tempz);
                }
            }
        }

        else
        {
            /*
               Things are a little more complicated in this case.  The
               variable is var(i,j), but i and j are themselves functions
               (eg. var(var(0,1),8)).  This warrants a little more care.

               Cautionary note: Local variables must be absolute - ie.
               a variable of the form var(0,i) *must* have i an integer -
               otherwise, we'll just ignore it here.  In fact, we are
               *only* going to go to extra effort for global variables
               (ie. i>0, j>0).  The reason for this is that when i<=0 or
               j<=0 we assume that this is just a marker (for derivatives
               for example), and thus *cannot*, by definition, be
               complicated.  If this was not done, things would get *really
               ugly*.
            */

            if ( ( naive == 0 ) && ( i > 0 ) && ( j > 0 ) )
            {
                tempz = equation->SubFunctions;

                /*
                   xglobal_subvar_an is the following:

                   esplit(var(-10,-10),1,esplit(var(-20,-20),1,1,var(var(-10,-10),var(-20,-20))),var(var(-10,-10),var(-20,-20)))
                */

                if ( ( resulta = copyEquation(xglobal_subvar_an) ) == NULL )
                {
                    return EERR_BADCOPY;
                }

                /*
                   var(-10,-10) = "i"
                */

                if ( ( q = locsubVar(resulta,(equation->SubFunctions)[0],-10,-10,0) ) )
                {
                    return q;
                }

                /*
                   var(-20,-20) = "j"
                */

                if ( ( q = locsubVar(resulta,(equation->SubFunctions)[1],-20,-20,0) ) )
                {
                    return q;
                }

                resultb = (resulta->SubFunctions)[2];

                ((resulta->SubFunctions)[1])->multiplier = ResultInteger(i);
                ((resultb->SubFunctions)[1])->multiplier = ResultInteger(j);

                ((resultb->SubFunctions)[2])->CombinerFnNumber = var_replace->CombinerFnNumber;
                ((resultb->SubFunctions)[2])->multiplier       = var_replace->multiplier;

                k = (FunctionTable[var_replace->CombinerFnNumber]).NumLocInput;

                if ( k > 0 )
                {
                    if ( ( ((resultb->SubFunctions)[2])->SubFunctions = (MathsNode **) malloc(k*sizeof(MathsNode *)) ) == NULL )
                    {
                        ((resultb->SubFunctions)[2])->SubFunctions = tempz;

                        return EERR_MALLOCFAIL;
                    }

                    for ( l = 0 ; l < k ; l++ )
                    {
                        if ( ( (((resultb->SubFunctions)[2])->SubFunctions)[l] = copyEquation((var_replace->SubFunctions)[l]) ) == NULL )
                        {
                            if ( l > 0 )
                            {
                                for ( k = 0 ; k < l ; k++ )
                                {
                                    delEqn((((resultb->SubFunctions)[2])->SubFunctions)[k]);
                                }
                            }

                            free(((resultb->SubFunctions)[2])->SubFunctions);

                            ((resultb->SubFunctions)[2])->SubFunctions = tempz;

                            return EERR_BADCOPY;
                        }
                    }
                }

                else
                {
                    ((resultb->SubFunctions)[2])->SubFunctions = NULL;
                }

                equation->CombinerFnNumber = resulta->CombinerFnNumber;
                equation->multiplier       = mulGenVar(equation->multiplier,resulta->multiplier);
                equation->SubFunctions     = resulta->SubFunctions;

                /*
                   We *do not* want to free the subfunctions of resulta
                   here - they are being used in equation.
                */

                free(resulta); /* remember: we're using the subfunctions */
                delEqn(tempz[0]);
                delEqn(tempz[1]);
                free(tempz);
            }
        }
    }

    return 0;
}

int locsubVarRow(MathsNode *equation, MathsNode *var_replace, long i, int naive)
{
    long k,l;
    int q;
    GenVar tempi;
    MathsNode *expressj;
    MathsNode *resulta;
    MathsNode **tempz;
    MathsNode johnny_howard_is_destroying_this_country;

    /*
       replacing var(i,~) with ~ will not work directly - so instead
       we replace var(i,~) with 0+~.  The 0 will disappear after
       simplification, so the result will be as if we had just replaced
       var(i,~) with ~.
    */

    tempi = var_replace->multiplier;
    var_replace->multiplier = ResultOne();

    if ( ( var_replace->CombinerFnNumber                      == VAR_FN_NUMBER  ) &&
         ( ((var_replace->SubFunctions)[0])->CombinerFnNumber == ZERO_FN_NUMBER ) &&
         ( ((var_replace->SubFunctions)[1])->CombinerFnNumber == ZERO_FN_NUMBER )    )
    {
        var_replace->multiplier = tempi;

        johnny_howard_is_destroying_this_country = *xsubvar_spec_case_rep;
        johnny_howard_is_destroying_this_country.multiplier = mulGenVar(johnny_howard_is_destroying_this_country.multiplier,tempi);

        /*
           A note: it may seem that I have made a stupid mistake here
           (pointer assignment).  In this case, however, it does not
           matter.  Think about it.
        */

        return locsubVarRow(equation,&johnny_howard_is_destroying_this_country,i,naive);
    }

    var_replace->multiplier = tempi;

    k = (FunctionTable[equation->CombinerFnNumber]).NumLocInput;

    if ( k > 0 )
    {
        for ( l = 0 ; l < k ; l++ )
        {
            if ( ( q = locsubVarRow((equation->SubFunctions)[l],var_replace,i,naive) ) )
            {
                return q;
            }
        }
    }

    if ( ( equation->CombinerFnNumber == VAR_FN_NUMBER ) ||
         ( equation->CombinerFnNumber == X_FN_NUMBER   ) ||
         ( equation->CombinerFnNumber == Y_FN_NUMBER   ) ||
         ( equation->CombinerFnNumber == Z_FN_NUMBER   )    )
    {
        /*
           OK, we've hit a variable of some sort.  Find out which one.
        */

        switch ( equation->CombinerFnNumber )
        {
            case X_FN_NUMBER:
            {
                tempi = ResultOne();

                if ( ( expressj = copyEquation(xglobal_one) ) == NULL )
                {
                    return EERR_BADCOPY;
                }

                equation->SubFunctions = NULL;

                break;
            }

            case Y_FN_NUMBER:
            {
                tempi = ResultOne();

                if ( ( expressj = copyEquation(xglobal_two) ) == NULL )
                {
                    return EERR_BADCOPY;
                }

                equation->SubFunctions = NULL;

                break;
            }

            case Z_FN_NUMBER:
            {
                tempi = ResultOne();

                if ( ( expressj = copyEquation(xglobal_three) ) == NULL )
                {
                    return EERR_BADCOPY;
                }

                equation->SubFunctions = NULL;

                break;
            }

            default:
            {
                if ( ( q = Gen_evaluateEqnNull(&tempi,(equation->SubFunctions)[0],1,0) ) )
                {
                    return q;
                }

                if ( ( expressj = copyEquation((equation->SubFunctions)[1]) ) == NULL )
                {
                    return EERR_BADCOPY;
                }

                break;
            }
        }

        if ( isFiniteInteger(tempi) )
        {
            /*
               Easiest case - it's a variable where the relevant arguments is
               an integers (ie. known) - in this case it is straightforward
               to replace it.
            */

            if ( getInteger(tempi) == i )
            {
                /*
                   Save the current variable for future use.
                */

                tempz = equation->SubFunctions;

                equation->CombinerFnNumber = var_replace->CombinerFnNumber;
                equation->multiplier       = mulGenVar(equation->multiplier,var_replace->multiplier);

                k = (FunctionTable[equation->CombinerFnNumber]).NumLocInput;

                if ( k > 0 )
                {
                    if ( ( equation->SubFunctions = (MathsNode **) malloc(k*sizeof(MathsNode *)) ) == NULL )
                    {
                        equation->SubFunctions = tempz;

                        return EERR_MALLOCFAIL;
                    }

                    for ( l = 0 ; l < k ; l++ )
                    {
                        if ( ( (equation->SubFunctions)[l] = copyEquation((var_replace->SubFunctions)[l]) ) == NULL )
                        {
                            if ( l > 0 )
                            {
                                for ( k = 0 ; k < l ; k++ )
                                {
                                    delEqn((equation->SubFunctions)[k]);
                                }
                            }

                            free(equation->SubFunctions);

                            equation->SubFunctions = tempz;

                            return EERR_BADCOPY;
                        }

                        if ( ( q = locsubVar((equation->SubFunctions)[l],expressj,0,0,0) ) )
                        {
                            return q;
                        }
                    }
                }

                else
                {
                    equation->SubFunctions = NULL;
                }

                /*
                   Delete old variable.
                */

                if ( tempz != NULL )
                {
                    delEqn(tempz[0]);
                    delEqn(tempz[1]);
                    free(tempz);
                }
            }
        }

        else
        {
            /*
               Things are a little more complicated in this case.  The
               variable is var(i,j), but i and j are themselves functions
               (eg. var(var(0,1),8)).  This warrants a little more care.

               Cautionary note: Local variables must be absolute - ie.
               a variable of the form var(0,i) *must* have i an integer -
               otherwise, we'll just ignore it here.  In fact, we are
               *only* going to go to extra effort for global variables
               (ie. i>0, j>0).  The reason for this is that when i<=0 or
               j<=0 we assume that this is just a marker (for derivatives
               for example), and thus *cannot*, by definition, be
               complicated.  If this was not done, things get *really ugly*.
            */

            if ( ( i > 0 ) && ( naive == 0 ) )
            {
                tempz = equation->SubFunctions;

                /*
                   xglobal_subvar_ar is the following:

                   esplit(var(-10,-10),1,1,var(var(-10,-10),var(-20,-20)))
                */

                if ( ( resulta = copyEquation(xglobal_subvar_ar) ) == NULL )
                {
                    return EERR_BADCOPY;
                }

                /*
                   var(-10,-10) = "i"
                */

                if ( ( q = locsubVar(resulta,(equation->SubFunctions)[0],-10,-10,0) ) )
                {
                    return q;
                }

                /*
                   var(-20,-20) = "j"
                */

                if ( ( q = locsubVar(resulta,(equation->SubFunctions)[1],-20,-20,0) ) )
                {
                    return q;
                }

                ((resulta->SubFunctions)[1])->multiplier = ResultInteger(i);

                ((resulta->SubFunctions)[2])->CombinerFnNumber = var_replace->CombinerFnNumber;
                ((resulta->SubFunctions)[2])->multiplier       = var_replace->multiplier;

                k = (FunctionTable[var_replace->CombinerFnNumber]).NumLocInput;

                if ( k > 0 )
                {
                    if ( ( ((resulta->SubFunctions)[2])->SubFunctions = (MathsNode **) malloc(k*sizeof(MathsNode *)) ) == NULL )
                    {
                        ((resulta->SubFunctions)[2])->SubFunctions = tempz;

                        return EERR_MALLOCFAIL;
                    }

                    for ( l = 0 ; l < k ; l++ )
                    {
                        if ( ( (((resulta->SubFunctions)[2])->SubFunctions)[l] = copyEquation((var_replace->SubFunctions)[l]) ) == NULL )
                        {
                            if ( l > 0 )
                            {
                                for ( k = 0 ; k < l ; k++ )
                                {
                                    delEqn((((resulta->SubFunctions)[2])->SubFunctions)[k]);
                                }
                            }

                            free(((resulta->SubFunctions)[2])->SubFunctions);

                            ((resulta->SubFunctions)[2])->SubFunctions = tempz;

                            return EERR_BADCOPY;
                        }

                        if ( ( q = locsubVar((((resulta->SubFunctions)[2])->SubFunctions)[l],expressj,0,0,0) ) )
                        {
                            return q;
                        }
                    }
                }

                else
                {
                    ((resulta->SubFunctions)[2])->SubFunctions = NULL;
                }

                equation->CombinerFnNumber = resulta->CombinerFnNumber;
                equation->multiplier       = mulGenVar(equation->multiplier,resulta->multiplier);
                equation->SubFunctions     = resulta->SubFunctions;

                /*
                   We *do* *not* want to free the subfunctions of resulta
                   here - they are being used in equation.
                */

                free(resulta); /* remember: we're using the subfunctions */
                delEqn(tempz[0]);
                delEqn(tempz[1]);
                free(tempz);
            }
        }

        delEqn(expressj);
    }

    return 0;
}

int locsubVarCol(MathsNode *equation, MathsNode *var_replace, long j, int naive)
{
    long k,l;
    int q;
    GenVar tempj;
    GenVar tempi;
    MathsNode *expressi;
    MathsNode *resulta;
    MathsNode **tempz;
    MathsNode johnny_howard_is_destroying_this_country;

    /*
       replacing var(i,~) with ~ will not work directly - so instead
       we replace var(i,~) with 0+~.  The 0 will disappear after
       simplification, so the result will be as if we had just replaced
       var(i,~) with ~.
    */

    tempi = var_replace->multiplier;
    var_replace->multiplier = ResultOne();

    if ( ( var_replace->CombinerFnNumber                      == VAR_FN_NUMBER  ) &&
         ( ((var_replace->SubFunctions)[0])->CombinerFnNumber == ZERO_FN_NUMBER ) &&
         ( ((var_replace->SubFunctions)[1])->CombinerFnNumber == ZERO_FN_NUMBER )    )
    {
        var_replace->multiplier = tempi;

        johnny_howard_is_destroying_this_country = *xsubvar_spec_case_rep;
        johnny_howard_is_destroying_this_country.multiplier = mulGenVar(johnny_howard_is_destroying_this_country.multiplier,tempi);

        /*
           A note: it may seem that I have made a stupid mistake here
           (pointer assignment).  In this case, however, it does not
           matter.  Think about it.
        */

        return locsubVarCol(equation,&johnny_howard_is_destroying_this_country,j,naive);
    }

    var_replace->multiplier = tempi;

    k = (FunctionTable[equation->CombinerFnNumber]).NumLocInput;

    if ( k > 0 )
    {
        for ( l = 0 ; l < k ; l++ )
        {
            if ( ( q = locsubVarCol((equation->SubFunctions)[l],var_replace,j,naive) ) )
            {
                return q;
            }
        }
    }

    if ( ( equation->CombinerFnNumber == VAR_FN_NUMBER ) ||
         ( equation->CombinerFnNumber == X_FN_NUMBER   ) ||
         ( equation->CombinerFnNumber == Y_FN_NUMBER   ) ||
         ( equation->CombinerFnNumber == Z_FN_NUMBER   )    )
    {
        /*
           OK, we've hit a variable of some sort.  Find out which one.
        */

        switch ( equation->CombinerFnNumber )
        {
            case X_FN_NUMBER:
            {
                if ( ( expressi = copyEquation(xglobal_one) ) == NULL )
                {
                    return EERR_BADCOPY;
                }

                tempj = ResultOne();

                equation->SubFunctions = NULL;

                break;
            }

            case Y_FN_NUMBER:
            {
                if ( ( expressi = copyEquation(xglobal_two) ) == NULL )
                {
                    return EERR_BADCOPY;
                }

                tempj = ResultOne();

                equation->SubFunctions = NULL;

                break;
            }

            case Z_FN_NUMBER:
            {
                if ( ( expressi = copyEquation(xglobal_three) ) == NULL )
                {
                    return EERR_BADCOPY;
                }

                tempj = ResultOne();

                equation->SubFunctions = NULL;

                break;
            }

            default:
            {
                if ( ( expressi = copyEquation((equation->SubFunctions)[0]) ) == NULL )
                {
                    return EERR_BADCOPY;
                }

                if ( ( q = Gen_evaluateEqnNull(&tempj,(equation->SubFunctions)[1],1,0) ) )
                {
                    return q;
                }

                break;
            }
        }

        if ( isFiniteInteger(tempj) )
        {
            /*
               Easiest case - it's a variable where the relevant arguments is
               an integers (ie. known) - in this case it is straightforward
               to replace it.
            */

            if ( getInteger(tempj) == j )
            {
                /*
                   Save the current variable for future use.
                */

                tempz = equation->SubFunctions;

                equation->CombinerFnNumber = var_replace->CombinerFnNumber;
                equation->multiplier       = mulGenVar(equation->multiplier,var_replace->multiplier);

                k = (FunctionTable[equation->CombinerFnNumber]).NumLocInput;

                if ( k > 0 )
                {
                    if ( ( equation->SubFunctions = (MathsNode **) malloc(k*sizeof(MathsNode *)) ) == NULL )
                    {
                        equation->SubFunctions = tempz;

                        return EERR_MALLOCFAIL;
                    }

                    for ( l = 0 ; l < k ; l++ )
                    {
                        if ( ( (equation->SubFunctions)[l] = copyEquation((var_replace->SubFunctions)[l]) ) == NULL )
                        {
                            if ( l > 0 )
                            {
                                for ( k = 0 ; k < l ; k++ )
                                {
                                    delEqn((equation->SubFunctions)[k]);
                                }
                            }

                            free(equation->SubFunctions);

                            equation->SubFunctions = tempz;

                            return EERR_BADCOPY;
                        }

                        if ( ( q = locsubVar((equation->SubFunctions)[l],expressi,0,0,0) ) )
                        {
                            return q;
                        }
                    }
                }

                else
                {
                    equation->SubFunctions = NULL;
                }

                /*
                   Delete old variable.
                */

                if ( tempz != NULL )
                {
                    delEqn(tempz[0]);
                    delEqn(tempz[1]);
                    free(tempz);
                }
            }
        }

        else
        {
            /*
               Things are a little more complicated in this case.  The
               variable is var(i,j), but i and j are themselves functions
               (eg. var(var(0,1),8)).  This warrants a little more care.

               Cautionary note: Local variables must be absolute - ie.
               a variable of the form var(0,i) *must* have i an integer -
               otherwise, we'll just ignore it here.  In fact, we are
               *only* going to go to extra effort for global variables
               (ie. i>0, j>0).  The reason for this is that when i<=0 or
               j<=0 we assume that this is just a marker (for derivatives
               for example), and thus *cannot*, by definition, be
               complicated.  If this was not done, things get *really ugly*.
            */

            if ( ( j > 0 ) && ( naive == 0 ) )
            {
                tempz = equation->SubFunctions;

                /*
                   xglobal_subvar_ac is the following:

                   esplit(var(-20,-20),1,1,var(var(-10,-10),var(-20,-20)))
                */

                if ( ( resulta = copyEquation(xglobal_subvar_ac) ) == NULL )
                {
                    return EERR_BADCOPY;
                }

                /*
                   var(-10,-10) = "i"
                */

                if ( ( q = locsubVar(resulta,(equation->SubFunctions)[0],-10,-10,0) ) )
                {
                    return q;
                }

                /*
                   var(-20,-20) = "j"
                */

                if ( ( q = locsubVar(resulta,(equation->SubFunctions)[1],-20,-20,0) ) )
                {
                    return q;
                }

                ((resulta->SubFunctions)[1])->multiplier = ResultInteger(j);

                ((resulta->SubFunctions)[2])->CombinerFnNumber = var_replace->CombinerFnNumber;
                ((resulta->SubFunctions)[2])->multiplier       = var_replace->multiplier;

                k = (FunctionTable[var_replace->CombinerFnNumber]).NumLocInput;

                if ( k > 0 )
                {
                    if ( ( ((resulta->SubFunctions)[2])->SubFunctions = (MathsNode **) malloc(k*sizeof(MathsNode *)) ) == NULL )
                    {
                        ((resulta->SubFunctions)[2])->SubFunctions = tempz;

                        return EERR_MALLOCFAIL;
                    }

                    for ( l = 0 ; l < k ; l++ )
                    {
                        if ( ( (((resulta->SubFunctions)[2])->SubFunctions)[l] = copyEquation((var_replace->SubFunctions)[l]) ) == NULL )
                        {
                            if ( l > 0 )
                            {
                                for ( k = 0 ; k < l ; k++ )
                                {
                                    delEqn((((resulta->SubFunctions)[2])->SubFunctions)[k]);
                                }
                            }

                            free(((resulta->SubFunctions)[2])->SubFunctions);

                            ((resulta->SubFunctions)[2])->SubFunctions = tempz;

                            return EERR_BADCOPY;
                        }

                        if ( ( q = locsubVar((((resulta->SubFunctions)[2])->SubFunctions)[l],expressi,0,0,0) ) )
                        {
                            return q;
                        }
                    }
                }

                else
                {
                    ((resulta->SubFunctions)[2])->SubFunctions = NULL;
                }

                equation->CombinerFnNumber = resulta->CombinerFnNumber;
                equation->multiplier       = mulGenVar(equation->multiplier,resulta->multiplier);
                equation->SubFunctions     = resulta->SubFunctions;

                /*
                   We *do* *not* want to free the subfunctions of resulta
                   here - they are being used in equation.
                */

                free(resulta); /* remember: we're using the subfunctions */
                delEqn(tempz[0]);
                delEqn(tempz[1]);
                free(tempz);
            }
        }

        delEqn(expressi);
    }

    return 0;
}




/*********************************************************************

   Functions: differentiation and integration
   ==========================================

*********************************************************************/

MathsNode *makeDeriv(MathsNode *start, long i, long j)
{
    MathsNode *result;

    if ( ( result = fast_makeDeriv(start,i,j) ) == NULL )
    {
        return NULL;
    }

    if ( simplifyEqn(result) )
    {
        delEqn(result);

        return NULL;
    }

    return result;
}

MathsNode *MatrixMakeDeriv(MathsNode *start, long i, long j, long k, long l)
{
    MathsNode *result;

    if ( ( result = fast_MatrixMakeDeriv(start,i,j,k,l) ) == NULL )
    {
        return NULL;
    }

    if ( simplifyEqn(result) )
    {
        delEqn(result);

        return NULL;
    }

    return result;
}

MathsNode *RowMakeDeriv(MathsNode *start, long i, long j, long k)
{
    MathsNode *result;

    if ( ( result = fast_RowMakeDeriv(start,i,j,k) ) == NULL )
    {
        return NULL;
    }

    if ( simplifyEqn(result) )
    {
        delEqn(result);

        return NULL;
    }

    return result;
}

MathsNode *ColMakeDeriv(MathsNode *start, long i, long j, long k)
{
    MathsNode *result;

    if ( ( result = fast_ColMakeDeriv(start,i,j,k) ) == NULL )
    {
        return NULL;
    }

    if ( simplifyEqn(result) )
    {
        delEqn(result);

        return NULL;
    }

    return result;
}

MathsNode *fast_makeDeriv(MathsNode *start, long i, long j)
{
    MathsNode *result;
    GenArg locGlobInput;
    GenVar dummy;

    if ( ( i <= 0 ) || ( j <= 0 ) )
    {
        return NULL;
    }

    locGlobInput.argContents           = NULL;
    locGlobInput.ArgEvaluationFunction = &ArgEvalNull;

    if ( alloc_genarg_locals(&locGlobInput) )
    {
        return NULL;
    }

    /*
       First thing we need to do is mark all the local variables that
       are used in this equation.  Having done this, the derivative
       function will be able to allocate unused local variables, where
       necessary, to derivatives.

       That is all this line does - it just evaluates the expression
       (with touch_all = 1, so all branches are evaluated) so that
       this sets the flags for local variables that are used in this
       expression.
    */

    locevaluateEqn(&dummy,start,locGlobInput,1,1,0);

    /*
       Having done that, we can proceed to find the derivative.
    */

    result = locmakeDeriv(start,i,j,locGlobInput);

    dealloc_genarg_locals(&locGlobInput);

    return result;
}

MathsNode *fast_MatrixMakeDeriv(MathsNode *start, long i, long j, long k, long l)
{
    MathsNode *result;
    GenArg locGlobInput;
    GenVar dummy;

    if ( ( i <= 0 ) || ( j <= 0 ) )
    {
        return NULL;
    }

    locGlobInput.argContents           = NULL;
    locGlobInput.ArgEvaluationFunction = &ArgEvalNull;

    if ( alloc_genarg_locals(&locGlobInput) )
    {
        return NULL;
    }

    /*
       First thing we need to do is mark all the local variables that
       are used in this equation.  Having done this, the derivative
       function will be able to allocate unused local variables, where
       necessary, to derivatives.

       That is all this line does - it just evaluates the expression
       (with touch_all = 1, so all branches are evaluated) so that
       this sets the flags for local variables that are used in this
       expression.
    */

    locevaluateEqn(&dummy,start,locGlobInput,1,1,0);

    /*
       Having done that, we can proceed to find the derivative.
    */

    result = locMatrixMakeDeriv(start,i,j,k,l,locGlobInput);

    dealloc_genarg_locals(&locGlobInput);

    return result;
}

MathsNode *fast_RowMakeDeriv(MathsNode *start, long i, long j, long k)
{
    MathsNode *result;
    GenArg locGlobInput;
    GenVar dummy;

    if ( ( i <= 0 ) || ( j <= 0 ) )
    {
        return NULL;
    }

    locGlobInput.argContents           = NULL;
    locGlobInput.ArgEvaluationFunction = &ArgEvalNull;

    if ( alloc_genarg_locals(&locGlobInput) )
    {
        return NULL;
    }

    /*
       First thing we need to do is mark all the local variables that
       are used in this equation.  Having done this, the derivative
       function will be able to allocate unused local variables, where
       necessary, to derivatives.

       That is all this line does - it just evaluates the expression
       (with touch_all = 1, so all branches are evaluated) so that
       this sets the flags for local variables that are used in this
       expression.
    */

    locevaluateEqn(&dummy,start,locGlobInput,1,1,0);

    /*
       Having done that, we can proceed to find the derivative.
    */

    result = locRowMakeDeriv(start,i,j,k,locGlobInput);

    dealloc_genarg_locals(&locGlobInput);

    return result;
}

MathsNode *fast_ColMakeDeriv(MathsNode *start, long i, long j, long k)
{
    MathsNode *result;
    GenArg locGlobInput;
    GenVar dummy;

    if ( ( i <= 0 ) || ( j <= 0 ) )
    {
        return NULL;
    }

    locGlobInput.argContents           = NULL;
    locGlobInput.ArgEvaluationFunction = &ArgEvalNull;

    if ( alloc_genarg_locals(&locGlobInput) )
    {
        return NULL;
    }

    /*
       First thing we need to do is mark all the local variables that
       are used in this equation.  Having done this, the derivative
       function will be able to allocate unused local variables, where
       necessary, to derivatives.

       That is all this line does - it just evaluates the expression
       (with touch_all = 1, so all branches are evaluated) so that
       this sets the flags for local variables that are used in this
       expression.
    */

    locevaluateEqn(&dummy,start,locGlobInput,1,1,0);

    /*
       Having done that, we can proceed to find the derivative.
    */

    result = locColMakeDeriv(start,i,j,k,locGlobInput);

    dealloc_genarg_locals(&locGlobInput);

    return result;
}

MathsNode *locmakeDeriv(MathsNode *start, long i, long j, GenArg GlobInput)
{
    MathsNode *result;
    MathsNode *temp;

    if ( ( result = locMatrixMakeDeriv(start,-50,-50,-60,-60,GlobInput) ) == NULL )
    {
        return NULL;
    }

    if ( ( temp = copyEquation(xglobal_one) ) == NULL )
    {
        delEqn(result);

        return NULL;
    }

    temp->multiplier = ResultInteger(i);

    if ( locsubVar(result,temp,-50,-50,0) )
    {
        return NULL;
    }

    temp->multiplier = ResultInteger(j);

    if ( locsubVar(result,temp,-60,-60,0) )
    {
        return NULL;
    }

    delEqn(temp);

    return result;
}

MathsNode *locRowMakeDeriv(MathsNode *start, long i, long j, long k, GenArg GlobInput)
{
    MathsNode *result;
    MathsNode *temp;

    if ( ( result = locMatrixMakeDeriv(start,-50,-50,j,k,GlobInput) ) == NULL )
    {
        return NULL;
    }

    if ( ( temp = copyEquation(xglobal_one) ) == NULL )
    {
        delEqn(result);

        return NULL;
    }

    temp->multiplier = ResultInteger(i);

    if ( locsubVar(result,temp,-50,-50,0) )
    {
        return NULL;
    }

    delEqn(temp);

    return result;
}

MathsNode *locColMakeDeriv(MathsNode *start, long i, long j, long k, GenArg GlobInput)
{
    MathsNode *result;
    MathsNode *temp;

    if ( ( result = locMatrixMakeDeriv(start,i,j,-60,-60,GlobInput) ) == NULL )
    {
        return NULL;
    }

    if ( ( temp = copyEquation(xglobal_one) ) == NULL )
    {
        delEqn(result);

        return NULL;
    }

    temp->multiplier = ResultInteger(k);

    if ( locsubVar(result,temp,-60,-60,0) )
    {
        return NULL;
    }

    delEqn(temp);

    return result;
}

MathsNode *locMatrixMakeDeriv(MathsNode *start, long i, long j, long k, long l, GenArg GlobInput)
{
    MathsNode *result;
    MathsNode *temp;
    MathsNode *tempb;
    long mmm,vvv;

    /*
       Take a copy of the appropriate derivative template.
    */

    if ( ( result = copyEquation((FunctionTable[start->CombinerFnNumber]).Deriv_equ) ) == NULL )
    {
        return NULL;
    }

    if ( ( temp = copyEquation(xglobal_one) ) == NULL )
    {
        return NULL;
    }

    /*
       var(-3,-mmm) are placeholders for local variables - replace as
       required.
    */

    if ( (FunctionTable[start->CombinerFnNumber]).NumLocDerivVar )
    {
        /*
           OK - there are local variables that need to be made concrete.
           Loop through these, find the first unused local, and allocate
           this to them.  After doing this, touch the local variable to
           make sure that it doesn't get used multiple times.
        */

        for ( mmm = 1 ; mmm <= (FunctionTable[start->CombinerFnNumber]).NumLocDerivVar ; mmm++ )
        {
            if ( ( vvv = find_first_unused_loc_var(GlobInput) ) <= 0 )
            {
                return NULL;
            }

            touch_loc_var(vvv,GlobInput);

            temp->multiplier = ResultInteger(vvv);

            if ( locsubVar(result,temp,-3,-mmm,0) )
            {
                return NULL;
            }
        }
    }

    delEqn(temp);

    /*
       var(-4,-1) = var(i,j)
       var(-4,-2) = var(k,l)
    */

    if ( ( temp = copyEquation(xvar_generic) ) == NULL )
    {
        return NULL;
    }

    ((temp->SubFunctions)[0])->multiplier = ResultInteger(i);
    ((temp->SubFunctions)[1])->multiplier = ResultInteger(j);

    if ( locsubVar(result,temp,-4,-1,0) )
    {
        return NULL;
    }

    ((temp->SubFunctions)[0])->multiplier = ResultInteger(k);
    ((temp->SubFunctions)[1])->multiplier = ResultInteger(l);

    if ( locsubVar(result,temp,-4,-2,0) )
    {
        return NULL;
    }

    delEqn(temp);

    /*
       Now substitute subfunctions and subfunction variables as needed.
    */

    if ( (FunctionTable[start->CombinerFnNumber]).NumLocInput )
    {
        for ( mmm = 1 ; mmm <= (FunctionTable[start->CombinerFnNumber]).NumLocInput ; mmm++ )
        {
            if ( locsubVar(result,(start->SubFunctions)[mmm-1],-1,-mmm,0) )
            {
                return NULL;
            }

            if ( ( tempb = locMatrixMakeDeriv((start->SubFunctions)[mmm-1],i,j,k,l,GlobInput) ) == NULL )
            {
                return NULL;
            }

            if ( locsubVar(result,tempb,-2,-mmm,0) )
            {
                return NULL;
            }

            delEqn(tempb);
        }
    }

    if ( result != NULL )
    {
        result->multiplier = mulGenVar(result->multiplier,start->multiplier);
    }

    return result;
}




/*********************************************************************

Functions: evaluation
=====================

*********************************************************************/

double evaluateEqnNull(MathsNode *equation)
{
    double result;

    evaluateEqnNull_e(&result,equation);

    return result;
}

int evaluateEqnNull_e(double *answer, MathsNode *equation)
{
    GenVar result;
    int errcode;

    if ( ( errcode = Gen_evaluateEqnNull(&result,equation,0,0) ) )
    {
        *answer = GSL_NAN;

        return errcode;
    }

    *answer = getReal(result);

    return 0;
}

int Gen_evaluateEqnNull(GenVar *result, MathsNode *equation, int report_change, int simplify_onthefly)
{
    GenArg locGlobInput;
    int q;

    locGlobInput.argContents           = NULL;
    locGlobInput.ArgEvaluationFunction = &ArgEvalNull;

    /*
       Allocate local variables (or at least try).
    */

    if ( ( q = alloc_genarg_locals(&locGlobInput) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(result,equation,locGlobInput,0,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    dealloc_genarg_locals(&locGlobInput);

    return 0;
}

int ArgEvalNull(GenVar *result, void *argContents, long i, long j)
{
    *result = ResultUncalced();

    return 0;

    /*
       stuff to keep gcc from whinging.
    */

    argContents = NULL;
    i = j;
}

double evaluateEqnMatrix(MathsNode *equation, double **GlobInput)
{
    double result;

    evaluateEqnMatrix_e(&result,equation,GlobInput);

    return result;
}

int evaluateEqnMatrix_e(double *answer, MathsNode *equation, double **GlobInput)
{
    GenVar result;
    int errcode;

    if ( ( errcode = Gen_evaluateEqnMatrix(&result,equation,GlobInput) ) )
    {
        *answer = GSL_NAN;

        return errcode;
    }

    *answer = getReal(result);

    return 0;
}

int Gen_evaluateEqnMatrix(GenVar *result, MathsNode *equation, double **GlobInput)
{
    GenArg locGlobInput;
    int q;

    locGlobInput.argContents           = (void *) GlobInput;
    locGlobInput.ArgEvaluationFunction = &ArgEvalDoubleDouble;

    if ( ( q = alloc_genarg_locals(&locGlobInput) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(result,equation,locGlobInput,0,0,0) ) )
    {
        return q;
    }

    dealloc_genarg_locals(&locGlobInput);

    return 0;
}

int ArgEvalDouble(GenVar *result, void *argContents, long i, long j)
{
    double *contents;

    contents = (double *) argContents;

    if ( ( i == 1 ) && ( contents != NULL ) )
    {
        *result = ResultReal(contents[j-1]);
    }

    else
    {
        return EERR_BADARGS;
    }

    return 0;
}

double evaluateEqnVector(MathsNode *equation, double *GlobInput)
{
    double result;

    evaluateEqnVector_e(&result,equation,GlobInput);

    return result;
}

int evaluateEqnVector_e(double *answer, MathsNode *equation, double *GlobInput)
{
    GenVar result;
    int errcode;

    if ( ( errcode = Gen_evaluateEqnVector(&result,equation,GlobInput) ) )
    {
        *answer = GSL_NAN;

        return errcode;
    }

    *answer = getReal(result);

    return 0;
}

int Gen_evaluateEqnVector(GenVar *result, MathsNode *equation, double *GlobInput)
{
    GenArg locGlobInput;
    int q;

    locGlobInput.argContents           = (void *) GlobInput;
    locGlobInput.ArgEvaluationFunction = &ArgEvalDouble;

    if ( ( q = alloc_genarg_locals(&locGlobInput) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(result,equation,locGlobInput,0,0,0) ) )
    {
        return q;
    }

    dealloc_genarg_locals(&locGlobInput);

    return 0;
}

int ArgEvalDoubleDouble(GenVar *result, void *argContents, long i, long j)
{
    double **contents;

    contents = (double **) argContents;

    if ( contents != NULL )
    {
        *result = ResultReal(contents[i-1][j-1]);
    }

    else
    {
        return EERR_BADARGS;
    }

    return 0;
}

double evaluateEqnList(MathsNode *equation, int numargs, ...)
{
    double result;
    int j;
    GenArgList *temp;
    GenArg locGlobInput;
    GenVar xxx;
    va_list ap;

    /*
       NB: b/c of variable arg length stuff, we cannot exit this function
           prematurely.
    */

    va_start(ap, numargs);

    if ( numargs == 0 )
    {
        /*
           There are no arguments, so this is the same as the Null case.
        */

        result = evaluateEqnNull(equation);
    }

    else
    {
        if ( ( temp = (GenArgList *) malloc(sizeof(GenArgList)) ) == NULL )
        {
            result = GSL_NAN;
        }

        else
        {
            if ( ( temp->argContents = (double *) malloc(numargs*sizeof(double)) ) == NULL )
            {
                result = GSL_NAN;
            }

            else
            {
                temp->size = numargs;
    
                for ( j = 1 ; j <= numargs ; j++ )
                {
                    (temp->argContents)[j-1] = va_arg(ap,double);
                }

                locGlobInput.argContents           = (void *) temp;
                locGlobInput.ArgEvaluationFunction = &ArgEvalList;

                if ( alloc_genarg_locals(&locGlobInput) )
                {
                    result = GSL_NAN;
                }

                else
                {
                    if ( locevaluateEqn(&xxx,equation,locGlobInput,0,0,0) )
                    {
                        result = GSL_NAN;
                    }

                    else
                    {
                        result = getReal(xxx);

                        dealloc_genarg_locals(&locGlobInput);

                        free(temp->argContents);
                        free(temp);
                    }
                }
            }
        }
    }

    va_end(ap);

    return result;
}

int evaluateEqnList_e(double *answer, MathsNode *equation, int numargs, ...)
{
    int j;
    int errcode;
    GenArgList *temp;
    GenArg locGlobInput;
    GenVar xxx;
    va_list ap;

    /*
       NB: b/c of variable arg length stuff, we cannot exit this function
           prematurely.
    */

    va_start(ap, numargs);

    if ( numargs == 0 )
    {
        /*
           There are no arguments, so this is the same as the Null case.
        */

        errcode = evaluateEqnNull_e(answer,equation);
    }

    else
    {
        if ( ( temp = (GenArgList *) malloc(sizeof(GenArgList)) ) == NULL )
        {
            *answer = GSL_NAN;
            errcode = EERR_MALLOCFAIL;
        }

        else
        {
            if ( ( temp->argContents = (double *) malloc(numargs*sizeof(double)) ) == NULL )
            {
                *answer = GSL_NAN;
                errcode = EERR_MALLOCFAIL;
            }

            else
            {
                temp->size = numargs;
    
                for ( j = 1 ; j <= numargs ; j++ )
                {
                    (temp->argContents)[j-1] = va_arg(ap,double);
                }

                locGlobInput.argContents           = (void *) temp;
                locGlobInput.ArgEvaluationFunction = &ArgEvalList;

                if ( alloc_genarg_locals(&locGlobInput) )
                {
                    *answer = GSL_NAN;
                    errcode = EERR_MALLOCFAIL;
                }

                else
                {
                    if ( ( errcode = locevaluateEqn(&xxx,equation,locGlobInput,0,0,0) ) )
                    {
                        *answer = GSL_NAN;
                    }

                    else
                    {
                        *answer = getReal(xxx);
                        errcode = 0;

                        dealloc_genarg_locals(&locGlobInput);

                        free(temp->argContents);
                        free(temp);
                    }
                }
            }
        }
    }

    va_end(ap);

    return errcode;
}

int ArgEvalList(GenVar *result, void *tempa, long i, long j)
{
    GenArgList *temp;

    temp = (GenArgList *) tempa;

    if ( ( i == 1 ) && ( j >= 1 ) && ( j <= temp->size ) )
    {
        *result = ResultReal((temp->argContents)[j-1]);
    }

    else
    {
        return EERR_BADARGS;
    }

    return 0;
}

double evaluateEqnGeneral(MathsNode *equation, double (*DoubleArgEvaluationFunction)(void *argContents, long i, long j), void *GlobInput)
{
    double result;

    evaluateEqnGeneral_e(&result,equation,DoubleArgEvaluationFunction,GlobInput);

    return result;
}

int evaluateEqnGeneral_e(double *answer, MathsNode *equation, double (*DoubleArgEvaluationFunction)(void *argContents, long i, long j), void *GlobInput)
{
    GenVar result;
    int errcode;

    if ( ( errcode = Gen_evaluateEqnGeneral(&result,equation,DoubleArgEvaluationFunction,GlobInput) ) )
    {
        *answer = GSL_NAN;

        return errcode;
    }

    *answer = getReal(result);

    return 0;
}

int Gen_evaluateEqnGeneral(GenVar *result, MathsNode *equation, double (*DoubleArgEvaluationFunction)(void *argContents, long i, long j), void *GlobInput)
{
    GenArg locGlobInput;
    GenArgGeneral *temp;
    int q;

    if ( ( temp = (GenArgGeneral *) malloc(sizeof(GenArgGeneral)) ) == NULL )
    {
        return EERR_MALLOCFAIL;
    }

    temp->argContents           = GlobInput;
    temp->ArgEvaluationFunction = DoubleArgEvaluationFunction;

    locGlobInput.argContents           = (void *) temp;
    locGlobInput.ArgEvaluationFunction = &ArgEvalGeneral;

    if ( ( q = alloc_genarg_locals(&locGlobInput) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(result,equation,locGlobInput,0,0,0) ) )
    {
        return q;
    }

    dealloc_genarg_locals(&locGlobInput);

    free(temp);

    return 0;
}

int ArgEvalGeneral(GenVar *result, void *argContents, long i, long j)
{
    GenArgGeneral *temp;

    temp = (GenArgGeneral *) argContents;

    *result = ResultReal((temp->ArgEvaluationFunction)((temp->argContents),i,j));

    return 0;
}

int locevaluateEqn(GenVar *result, MathsNode *equation, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb;
    int q;

    /*
       NB: when simplification is occuring, the multiplier may change during
           the process of evaluating aa.  Therefore it would be inadvisable
           to glop this expression onto one line.
    */

    if ( ( q = ((FunctionTable[equation->CombinerFnNumber]).CalcOperation)(&aa,equation,GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    bb = equation->multiplier;

    *result = mulGenVar(aa,bb);

    return 0;
}

void make_number(MathsNode *equation, GenVar *result)
{
    int numsubs,i;

    /*
       We know that the equation has evaluated to the number certain result.
       (ie. result is not DATA_IS_UNCALCED).  What we need to do is work out
       a suitable substitute and, after deleting the old equation, replace
       the old result with the new.
    */

    numsubs = (FunctionTable[equation->CombinerFnNumber]).NumLocInput;

    if ( numsubs > 0 )
    {
        for ( i = 0 ; i < numsubs ; i++ )
        {
            delEqn((equation->SubFunctions)[i]);
        }

        free(equation->SubFunctions);
    }

    *result = mulGenVar(*result,equation->multiplier);

    switch ( result->DataType )
    {
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        case DATA_IS_ONE:
        {
            equation->CombinerFnNumber = ONE_FN_NUMBER;
            equation->multiplier       = *result;
            equation->SubFunctions     = NULL;

            *result = ResultOne();

            break;
        }

        case DATA_IS_ZERO:
        {
            equation->CombinerFnNumber = ZERO_FN_NUMBER;
            equation->multiplier       = ResultOne();
            equation->SubFunctions     = NULL;

            *result = ResultZero();

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            equation->CombinerFnNumber = PINF_FN_NUMBER;
            equation->multiplier       = ResultOne();
            equation->SubFunctions     = NULL;

            *result = ResultPosInfnty();

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            equation->CombinerFnNumber = NINF_FN_NUMBER;
            equation->multiplier       = ResultOne();
            equation->SubFunctions     = NULL;

            *result = ResultNegInfnty();

            break;
        }

        case DATA_IS_FIN_INDET:
        {
            equation->CombinerFnNumber = FIND_FN_NUMBER;
            equation->multiplier       = ResultOne();
            equation->SubFunctions     = NULL;

            *result = ResultFinIndet();

            break;
        }

        case DATA_IS_INF_INDET:
        {
            equation->CombinerFnNumber = IIND_FN_NUMBER;
            equation->multiplier       = ResultOne();
            equation->SubFunctions     = NULL;

            *result = ResultInfIndet();

            break;
        }

        case DATA_IS_ERROR:
        default:
        {
            equation->CombinerFnNumber = ERR_FN_NUMBER;
            equation->multiplier       = ResultOne();
            equation->SubFunctions     = (MathsNode **) malloc(sizeof(MathsNode *));

            *result = ResultError(GVAR_EUNKNOWN);

            break;
        }
    }

    return;
}




/*********************************************************************

   Functions: hacking
   ==================

*********************************************************************/

int initMaths(void)
{
    long i,j;
    char tempa[100];
    char *tempb;

    if ( ( xglobal_zero  = fast_parseMaths("0") ) == NULL ) { return 1; }
    if ( ( xglobal_one   = fast_parseMaths("1") ) == NULL ) { return 1; }
    if ( ( xglobal_two   = fast_parseMaths("2") ) == NULL ) { return 1; }
    if ( ( xglobal_three = fast_parseMaths("3") ) == NULL ) { return 1; }

    if ( ( tempb = deparseMaths_cstruct(xglobal_zero,"xglobal_zero") ) == NULL ) { return 1; }
    printf("%s\n",tempb);
    free(tempb);
    if ( ( tempb = deparseMaths_cstruct(xglobal_one,"xglobal_one") ) == NULL ) { return 1; }
    printf("%s\n",tempb);
    free(tempb);
    if ( ( tempb = deparseMaths_cstruct(xglobal_two,"xglobal_two") ) == NULL ) { return 1; }
    printf("%s\n",tempb);
    free(tempb);
    if ( ( tempb = deparseMaths_cstruct(xglobal_three,"xglobal_three") ) == NULL ) { return 1; }
    printf("%s\n",tempb);
    free(tempb);

    if ( ( xglobal_subvar_an = fast_parseMaths("esplit(var(-10,-10),1,esplit(var(-20,-20),1,1,var(var(-10,-10),var(-20,-20))),var(var(-10,-10),var(-20,-20)))")                       ) == NULL ) { return 1; }
    if ( ( xglobal_subvar_ar = fast_parseMaths("esplit(var(-10,-10),1,1,var(var(-10,-10),var(-20,-20)))")                                                                             ) == NULL ) { return 1; }
    if ( ( xglobal_subvar_ac = fast_parseMaths("esplit(var(-20,-20),1,1,var(var(-10,-10),var(-20,-20)))")                                                                             ) == NULL ) { return 1; }
    if ( ( xglobal_subvar_b  = fast_parseMaths("esplit(var(-10,-10),var(-20,-20),esplit(var(-30,-30),var(-40,-40),1,var(var(-10,-10),var(-30,-30))),var(var(-10,-10),var(-30,-30)))") ) == NULL ) { return 1; }

    if ( ( tempb = deparseMaths_cstruct(xglobal_subvar_an,"xglobal_subvar_an") ) == NULL ) { return 1; }
    printf("%s\n",tempb);
    free(tempb);
    if ( ( tempb = deparseMaths_cstruct(xglobal_subvar_ar,"xglobal_subvar_ar") ) == NULL ) { return 1; }
    printf("%s\n",tempb);
    free(tempb);
    if ( ( tempb = deparseMaths_cstruct(xglobal_subvar_ac,"xglobal_subvar_ac") ) == NULL ) { return 1; }
    printf("%s\n",tempb);
    free(tempb);
    if ( ( tempb = deparseMaths_cstruct(xglobal_subvar_b,"xglobal_subvar_b") ) == NULL ) { return 1; }
    printf("%s\n",tempb);
    free(tempb);

    if ( ( xsubvar_spec_case_det = fast_parseMaths("var(0,0)")   ) == NULL ) { return 1; }
    if ( ( xsubvar_spec_case_rep = fast_parseMaths("0+var(0,0)") ) == NULL ) { return 1; }

    if ( ( tempb = deparseMaths_cstruct(xsubvar_spec_case_det,"xsubvar_spec_case_det") ) == NULL ) { return 1; }
    printf("%s\n",tempb);
    free(tempb);
    if ( ( tempb = deparseMaths_cstruct(xsubvar_spec_case_rep,"xsubvar_spec_case_rep") ) == NULL ) { return 1; }
    printf("%s\n",tempb);
    free(tempb);

    if ( ( xinteg_simplifier   = fast_parseMaths("var(-1,-1)*(var(-2,-2)-var(-3,-3))")     ) == NULL ) { return 1; }
    if ( ( xsumprod_simplifier = fast_parseMaths("eval(var(-1,-1),var(-2,-2),var(-3,-3))") ) == NULL ) { return 1; }

    if ( ( tempb = deparseMaths_cstruct(xinteg_simplifier,"xinteg_simplifier") ) == NULL ) { return 1; }
    printf("%s\n",tempb);
    free(tempb);
    if ( ( tempb = deparseMaths_cstruct(xsumprod_simplifier,"xsumprod_simplifier") ) == NULL ) { return 1; }
    printf("%s\n",tempb);
    free(tempb);

    if ( ( xsum_recurser  = fast_parseMaths("eval(var(-1,-1),var(-2,-2),var(-5,-5))+sum(var(-1,-1),var(-2,-2)+var(-4,-4),var(-3,-3),var(-4,-4),var(-5,-5))")  ) == NULL ) { return 1; }
    if ( ( xprod_recurser = fast_parseMaths("eval(var(-1,-1),var(-2,-2),var(-5,-5))*prod(var(-1,-1),var(-2,-2)+var(-4,-4),var(-3,-3),var(-4,-4),var(-5,-5))") ) == NULL ) { return 1; }

    if ( ( tempb = deparseMaths_cstruct(xsum_recurser,"xsum_recurser") ) == NULL ) { return 1; }
    printf("%s\n",tempb);
    free(tempb);
    if ( ( tempb = deparseMaths_cstruct(xprod_recurser,"xprod_recurser") ) == NULL ) { return 1; }
    printf("%s\n",tempb);
    free(tempb);

    if ( ( xvar_generic = fast_parseMaths("var(5,5)") ) == NULL ) { return 1; }

    if ( ( tempb = deparseMaths_cstruct(xvar_generic,"xvar_generic") ) == NULL ) { return 1; }
    printf("%s\n",tempb);
    free(tempb);

    for ( i = 1 ; i <= NUM_FUNCTIONS ; i++ )
    {
        if ( ( (FunctionTable[i-1]).Deriv_equ = fast_parseMaths((FunctionTable[i-1]).Deriv_string) ) == NULL ) { return 1; }

        sprintf(tempa,"%s_wderiv",(FunctionTable[i-1]).FnName);

        if ( ( tempb = deparseMaths_cstruct((FunctionTable[i-1]).Deriv_equ,tempa) ) == NULL ) { return 1; }

        for ( j = strlen(tempb)-1 ; j >= 0 ; j-- )
        {
            if ( tempb[j-1] == '\n' )
            {
                goto exit_point;
            }

            if ( tempb[j-1] == '=' )
            {
                tempb[j-1] = ' ';
            }
        }

        exit_point:

        tempb[j+0]  = '#'; /* M */
        tempb[j+1]  = 'd'; /* a */
        tempb[j+2]  = 'e'; /* t */
        tempb[j+3]  = 'f'; /* h */
        tempb[j+4]  = 'i'; /* s */
        tempb[j+5]  = 'n'; /* N */
        tempb[j+6]  = 'e'; /* o */
        tempb[j+7]  = ' '; /* d */
        tempb[j+8]  = ' '; /* e */
        tempb[j+9]  = ' '; /*   */
        tempb[j+10] = ' '; /* * */

        tempb[strlen(tempb)-2] = tempb[strlen(tempb)-1];

        printf("%s",tempb);

        free(tempb);
    }

    return 0;
}




/*********************************************************************

Local Variables
===============

*********************************************************************/

int alloc_genarg_locals(GenArg *context)
{
    long j;

    if ( ( context->loc_vars = (GenVar **) malloc(sizeof(GenVar *)) ) == NULL )
    {
        return EERR_MALLOCFAIL;
    }

    if ( ( context->loc_var_used_flags = (int **) malloc(sizeof(int *)) ) == NULL )
    {
        return EERR_MALLOCFAIL;
    }

    if ( ( (*(context->loc_vars)) = (GenVar *) malloc(MAX_NUM_LOC_VARS*sizeof(GenVar)) ) == NULL )
    {
        return EERR_MALLOCFAIL;
    }

    if ( ( (*(context->loc_var_used_flags)) = (int *) malloc(MAX_NUM_LOC_VARS*sizeof(int)) ) == NULL )
    {
        return EERR_MALLOCFAIL;
    }

    for ( j = 1 ; j <= MAX_NUM_LOC_VARS ; j++ )
    {
        (*(context->loc_var_used_flags))[j-1] = 0;
        (*(context->loc_vars))[j-1] = ResultUncalced();
    }

    return 0;
}

void dealloc_genarg_locals(GenArg *context)
{
    free((context->loc_vars)[0]);
    free((context->loc_var_used_flags)[0]);

    free(context->loc_vars);
    free(context->loc_var_used_flags);

    return;
}

long find_first_unused_loc_var(GenArg context)
{
    long i;

    for ( i = 1 ; i <= MAX_NUM_LOC_VARS ; i++ )
    {
        if ( (*(context.loc_var_used_flags))[i-1] == 0 )
        {
            return i;
        }
    }

    return 0;
}

int set_loc_var(long i, GenVar what, GenArg context)
{
    if ( ( i >= 1 ) && ( i <= MAX_NUM_LOC_VARS ) )
    {
        (*(context.loc_vars))[i-1] = what;

        (*(context.loc_var_used_flags))[i-1] = 1;

        return 0;
    }

    return -1;
}

int restore_loc_var(long i, GenVar what, GenArg context)
{
    if ( ( i >= 1 ) && ( i <= MAX_NUM_LOC_VARS ) )
    {
        (*(context.loc_vars))[i-1] = what;

        return 0;
    }

    return -1;
}

GenVar get_loc_var(long i, GenArg context)
{
    GenVar result;

    if ( ( i >= 1 ) && ( i <= MAX_NUM_LOC_VARS ) )
    {
        if ( (*(context.loc_var_used_flags))[i-1] )
        {
            result = (*(context.loc_vars))[i-1];

            return result;
        }

        else
        {
            result = ResultUncalced();

            return result;
        }
    }

    result = ResultError(GVAR_EBADLOCAL);

    return result;
}

int touch_loc_var(long i, GenArg context)
{
    if ( ( i >= 1 ) && ( i <= MAX_NUM_LOC_VARS ) )
    {
        (*(context.loc_var_used_flags))[i-1] = 1;

        return 0;
    }

    return -1;
}



























/*********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
*****************************      ***********************************
************************                ******************************
*********************                      ***************************
********************    Function section    **************************
********************    ================    **************************
*********************                      ***************************
************************                ******************************
*****************************      ***********************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
**********************************************************************
*********************************************************************/




#define REMOVE_GCC_WARNING_CALC         \
GlobInput.argContents = NULL;           \
touch_all = report_change;

int first_rand_flag = 1;

#define RAND_INIT_SEQUENCE                                              \
{                                                                       \
    if ( first_rand_flag )                                              \
    {                                                                   \
        first_rand_flag = 0;                                            \
                                                                        \
        srand(time(0));                                                 \
        RandomInitialise((rand()%31329),(rand()%30082));                \
    }                                                                   \
}

/*
HACK REMOVALS: simplification of sum, prod, integ removed as it was buggy

UNIMPLEMENTED SIMPLIFICATIONS:

test sum(m,i,j,k,x)   - if i == j then replace with trivially expanded x
test prod(m,i,j,k,x)  - if i == j then replace with trivially expanded x
test integ(m,a,b,c,x) - if a == b then replace with 0
                      - if x is a constant, simplify
 eval(m,i,x)      - always simplify trivially by substitution
*/



int Calc_one(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    (*result) = ResultOne();

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_zero(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    (*result) = ResultZero();

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_inf(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    (*result) = ResultPosInfnty();

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_ninf(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    (*result) = ResultNegInfnty();

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_FiniteIndet(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    (*result) = ResultFinIndet();

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_InfintIndet(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    (*result) = ResultInfIndet();

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_isError(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    (*result) = ResultError(GVAR_EUNKNOWN);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_var(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar ii,jj;
    long i,j,q;

    if ( ( q = locevaluateEqn(&ii,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&jj,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( isFiniteInteger(ii) && isFiniteInteger(jj) )
    {
        i = getInteger(ii);
        j = getInteger(jj);

        if ( i == 0 )
        {
            if ( j > 0 )
            {
                if ( simplify_onthefly )
                {
                    (*result) = ResultUncalced();
                }

                else
                {
                    (*result) = get_loc_var(j,GlobInput);

                    if ( touch_all )
                    {
                        touch_loc_var(j,GlobInput);
                    }
                }
            }

            else
            {
                if ( j < 0 )
                {
                    return EERR_INVALIDVAR;
                }

                else
                {
                    (*result) = ResultUncalced();
                }
            }
        }

        else
        {
            if ( ( i > 0 ) && ( j > 0 ) )
            {
                if ( ( q = (GlobInput.ArgEvaluationFunction)(result,(GlobInput.argContents),i,j) ) )
                {
                    return q;
                }
            }

            else
            {
                if ( ( i < 0 ) || ( j < 0 ) )
                {
                    return EERR_INVALIDVAR;
                }

                else
                {
                    (*result) = ResultUncalced();
                }
            }
        }
    }

    else
    {
        if ( ( ii.DataType == DATA_IS_UNCALCED ) ||
             ( jj.DataType == DATA_IS_UNCALCED )    )
        {
            (*result) = ResultUncalced();
        }

        else
        {
            return EERR_INVALIDVAR;
        }
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_x(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    int q;

    if ( ( q = (GlobInput.ArgEvaluationFunction)(result,(GlobInput.argContents),1,1) ) )
    {
        return q;
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_y(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    int q;

    if ( ( q = (GlobInput.ArgEvaluationFunction)(result,(GlobInput.argContents),1,2) ) )
    {
        return q;
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_z(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    int q;

    if ( ( q = (GlobInput.ArgEvaluationFunction)(result,(GlobInput.argContents),1,3) ) )
    {
        return q;
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_kronDelta(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar ii,jj;
    int q;

    if ( ( q = locevaluateEqn(&ii,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&jj,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_kronDelta(ii,jj);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_diracDelta(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar xx,yy;
    int q;

    if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_diracDelta(xx,yy);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_finDiracDel(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar xx,yy;
    int q;

    if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_finDiracDelta(xx,yy);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_split(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb,xx,yy;
    MathsNode **tempb;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&bb,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( isReal(aa) && isReal(bb) )
    {
        if ( isalessequalb(aa,bb) )
        {
            if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[3],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            (*result) = xx;
        }

        else
        {
            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[3],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            (*result) = yy;
        }
    }

    else
    {
        if ( ( aa.DataType == DATA_IS_UNCALCED ) ||
             ( bb.DataType == DATA_IS_UNCALCED )    )
        {
            if ( touch_all || simplify_onthefly )
            {
                if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }

                if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[3],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            (*result) = ResultUncalced();
        }

        else
        {
            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }

                if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[3],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            (*result) = ResultError(GVAR_GSL_EDOM);
        }
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }

        else
        {
            if ( isReal(aa) && isReal(bb) )
            {
                if ( isalessequalb(aa,bb) )
                {
                    tempb = _this->SubFunctions;

                    (tempb[2])->multiplier = mulGenVar(_this->multiplier,(tempb[2])->multiplier);

                    *_this = *(tempb[2]);

                    delEqn(tempb[0]);
                    delEqn(tempb[1]);
                    free(tempb[2]);
                    delEqn(tempb[3]);

                    free(tempb);
                }

                else
                {
                    tempb = _this->SubFunctions;

                    (tempb[3])->multiplier = mulGenVar(_this->multiplier,(tempb[3])->multiplier);

                    *_this = *(tempb[3]);

                    delEqn(tempb[0]);
                    delEqn(tempb[1]);
                    delEqn(tempb[2]);
                    free(tempb[3]);

                    free(tempb);
                }
            }

            else
            {
                if ( isReal(xx) && isReal(yy) )
                {
                    if ( isaequalb(xx,yy) )
                    {
                        make_number(_this,result);
                    }
                }
            }
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_esplit(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb,xx,yy;
    MathsNode **tempb;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&bb,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( isReal(aa) && isReal(bb) )
    {
        if ( isaequalb(aa,bb) )
        {
            if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[3],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            (*result) = xx;
        }

        else
        {
            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[3],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            (*result) = yy;
        }
    }

    else
    {
        if ( ( aa.DataType == DATA_IS_UNCALCED ) ||
             ( bb.DataType == DATA_IS_UNCALCED )    )
        {
            if ( touch_all || simplify_onthefly )
            {
                if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }

                if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[3],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            (*result) = ResultUncalced();
        }

        else
        {
            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }

                if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[3],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            (*result) = ResultError(GVAR_GSL_EDOM);
        }
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }

        else
        {
            if ( isReal(aa) && isReal(bb) )
            {
                if ( isaequalb(aa,bb) )
                {
                    tempb = _this->SubFunctions;

                    (tempb[2])->multiplier = mulGenVar(_this->multiplier,(tempb[2])->multiplier);

                    *_this = *(tempb[2]);

                    delEqn(tempb[0]);
                    delEqn(tempb[1]);
                    free(tempb[2]);
                    delEqn(tempb[3]);

                    free(tempb);
                }

                else
                {
                    tempb = _this->SubFunctions;

                    (tempb[3])->multiplier = mulGenVar(_this->multiplier,(tempb[3])->multiplier);

                    *_this = *(tempb[3]);

                    delEqn(tempb[0]);
                    delEqn(tempb[1]);
                    delEqn(tempb[2]);
                    free(tempb[3]);

                    free(tempb);
                }
            }

            else
            {
                if ( isReal(xx) && isReal(yy) )
                {
                    if ( isaequalb(xx,yy) )
                    {
                        make_number(_this,result);
                    }
                }
            }
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_isint(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,xx,yy;
    MathsNode **temp;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( aa.DataType == DATA_IS_UNCALCED )
    {
        if ( touch_all || simplify_onthefly )
        {
            if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }
        }

        (*result) = ResultUncalced();
    }

    else
    {
        if ( isInteger(aa) )
        {
            if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            (*result) = xx;
        }

        else
        {
            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            (*result) = yy;
        }
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }

        else
        {
            if ( aa.DataType != DATA_IS_UNCALCED )
            {
                if ( isInteger(aa) )
                {
                    temp = _this->SubFunctions;

                    (temp[1])->multiplier = mulGenVar(_this->multiplier,(temp[1])->multiplier);

                    *_this = *(temp[1]);

                    delEqn(temp[0]);
                    free(temp[1]);
                    delEqn(temp[2]);

                    free(temp);
                }

                else
                {
                    temp = _this->SubFunctions;

                    (temp[2])->multiplier = mulGenVar(_this->multiplier,(temp[2])->multiplier);

                    *_this = *(temp[2]);

                    delEqn(temp[0]);
                    delEqn(temp[1]);
                    free(temp[2]);

                    free(temp);
                }
            }

            else
            {
                if ( isReal(xx) && isReal(yy) )
                {
                    if ( isaequalb(xx,yy) )
                    {
                        *result = xx;

                        make_number(_this,result);
                    }
                }
            }
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_isreal(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,xx,yy;
    MathsNode **temp;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( aa.DataType == DATA_IS_UNCALCED )
    {
        if ( touch_all || simplify_onthefly )
        {
            if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }
        }

        (*result) = ResultUncalced();
    }

    else
    {
        if ( !isInteger(aa) && isReal(aa) )
        {
            if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            (*result) = xx;
        }

        else
        {
            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            (*result) = yy;
        }
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }

        else
        {
            if ( aa.DataType != DATA_IS_UNCALCED )
            {
                if ( !isInteger(aa) && isReal(aa) )
                {
                    temp = _this->SubFunctions;

                    (temp[1])->multiplier = mulGenVar(_this->multiplier,(temp[1])->multiplier);

                    *_this = *(temp[1]);

                    delEqn(temp[0]);
                    free(temp[1]);
                    delEqn(temp[2]);

                    free(temp);
                }

                else
                {
                    temp = _this->SubFunctions;

                    (temp[2])->multiplier = mulGenVar(_this->multiplier,(temp[2])->multiplier);

                    *_this = *(temp[2]);

                    delEqn(temp[0]);
                    delEqn(temp[1]);
                    free(temp[2]);

                    free(temp);
                }
            }

            else
            {
                if ( isReal(xx) && isReal(yy) )
                {
                    if ( isaequalb(xx,yy) )
                    {
                        (*result) = xx;

                        make_number(_this,result);
                    }
                }
            }
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_isanum(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,xx,yy;
    MathsNode **temp;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( aa.DataType == DATA_IS_UNCALCED )
    {
        if ( touch_all || simplify_onthefly )
        {
            if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }
        }

        (*result) = ResultUncalced();
    }

    else
    {
        if ( isReal(aa) )
        {
            if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            (*result) = xx;
        }

        else
        {
            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            (*result) = yy;
        }
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }

        else
        {
            if ( aa.DataType != DATA_IS_UNCALCED )
            {
                if ( isReal(aa) )
                {
                    temp = _this->SubFunctions;

                    (temp[1])->multiplier = mulGenVar(_this->multiplier,(temp[1])->multiplier);

                    *_this = *(temp[1]);

                    delEqn(temp[0]);
                    free(temp[1]);
                    delEqn(temp[2]);

                    free(temp);
                }

                else
                {
                    temp = _this->SubFunctions;

                    (temp[2])->multiplier = mulGenVar(_this->multiplier,(temp[2])->multiplier);

                    *_this = *(temp[2]);

                    delEqn(temp[0]);
                    delEqn(temp[1]);
                    free(temp[2]);

                    free(temp);
                }
            }

            else
            {
                if ( isReal(xx) && isReal(yy) )
                {
                    if ( isaequalb(xx,yy) )
                    {
                        (*result) = xx;

                        make_number(_this,result);
                    }
                }
            }
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_isfinint(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,xx,yy;
    MathsNode **temp;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( aa.DataType == DATA_IS_UNCALCED )
    {
        if ( touch_all || simplify_onthefly )
        {
            if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }
        }

        (*result) = ResultUncalced();
    }

    else
    {
        if ( isFiniteInteger(aa) )
        {
            if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            (*result) = xx;
        }

        else
        {
            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            (*result) = yy;
        }
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }

        else
        {
            if ( aa.DataType != DATA_IS_UNCALCED )
            {
                if ( isFiniteInteger(aa) )
                {
                    temp = _this->SubFunctions;

                    (temp[1])->multiplier = mulGenVar(_this->multiplier,(temp[1])->multiplier);

                    *_this = *(temp[1]);

                    delEqn(temp[0]);
                    free(temp[1]);
                    delEqn(temp[2]);

                    free(temp);
                }

                else
                {
                    temp = _this->SubFunctions;

                    (temp[2])->multiplier = mulGenVar(_this->multiplier,(temp[2])->multiplier);

                    *_this = *(temp[2]);

                    delEqn(temp[0]);
                    delEqn(temp[1]);
                    free(temp[2]);

                    free(temp);
                }
            }

            else
            {
                if ( isReal(xx) && isReal(yy) )
                {
                    if ( isaequalb(xx,yy) )
                    {
                        (*result) = xx;

                        make_number(_this,result);
                    }
                }
            }
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_isfinreal(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,xx,yy;
    MathsNode **temp;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( aa.DataType == DATA_IS_UNCALCED )
    {
        if ( touch_all || simplify_onthefly )
        {
            if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }
        }

        (*result) = ResultUncalced();
    }

    else
    {
        if ( !isFiniteInteger(aa) && isFiniteReal(aa) )
        {
            if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            (*result) = xx;
        }

        else
        {
            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            (*result) = yy;
        }
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }

        else
        {
            if ( aa.DataType != DATA_IS_UNCALCED )
            {
                if ( !isFiniteInteger(aa) && isFiniteReal(aa) )
                {
                    temp = _this->SubFunctions;

                    (temp[1])->multiplier = mulGenVar(_this->multiplier,(temp[1])->multiplier);

                    *_this = *(temp[1]);

                    delEqn(temp[0]);
                    free(temp[1]);
                    delEqn(temp[2]);

                    free(temp);
                }

                else
                {
                    temp = _this->SubFunctions;

                    (temp[2])->multiplier = mulGenVar(_this->multiplier,(temp[2])->multiplier);

                    *_this = *(temp[2]);

                    delEqn(temp[0]);
                    delEqn(temp[1]);
                    free(temp[2]);

                    free(temp);
                }
            }

            else
            {
                if ( isReal(xx) && isReal(yy) )
                {
                    if ( isaequalb(xx,yy) )
                    {
                        (*result) = xx;

                        make_number(_this,result);
                    }
                }
            }
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_isfinanum(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,xx,yy;
    MathsNode **temp;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( aa.DataType == DATA_IS_UNCALCED )
    {
        if ( touch_all || simplify_onthefly )
        {
            if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }
        }

        (*result) = ResultUncalced();
    }

    else
    {
        if ( isFiniteReal(aa) )
        {
            if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            (*result) = xx;
        }

        else
        {
            if ( touch_all )
            {
                if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            (*result) = yy;
        }
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }

        else
        {
            if ( aa.DataType != DATA_IS_UNCALCED )
            {
                if ( isFiniteReal(aa) )
                {
                    temp = _this->SubFunctions;

                    (temp[1])->multiplier = mulGenVar(_this->multiplier,(temp[1])->multiplier);

                    *_this = *(temp[1]);

                    delEqn(temp[0]);
                    free(temp[1]);
                    delEqn(temp[2]);

                    free(temp);
                }

                else
                {
                    temp = _this->SubFunctions;

                    (temp[2])->multiplier = mulGenVar(_this->multiplier,(temp[2])->multiplier);

                    *_this = *(temp[2]);

                    delEqn(temp[0]);
                    delEqn(temp[1]);
                    free(temp[2]);

                    free(temp);
                }
            }

            else
            {
                if ( isReal(xx) && isReal(yy) )
                {
                    if ( isaequalb(xx,yy) )
                    {
                        (*result) = xx;

                        make_number(_this,result);
                    }
                }
            }
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_add(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb;
    MathsNode **tempb;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&bb,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = addGenVar(aa,bb);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }

        else
        {
            if ( aa.DataType == DATA_IS_ZERO )
            {
                tempb = _this->SubFunctions;

                (tempb[1])->multiplier = mulGenVar(_this->multiplier,(tempb[1])->multiplier);

                *_this = *(tempb[1]);

                delEqn(tempb[0]);
                free(tempb[1]);

                free(tempb);
            }

            else
            {
                if ( bb.DataType == DATA_IS_ZERO )
                {
                    tempb = _this->SubFunctions;

                    (tempb[0])->multiplier = mulGenVar(_this->multiplier,(tempb[0])->multiplier);

                    *_this = *(tempb[0]);

                    free(tempb[0]);
                    delEqn(tempb[1]);

                    free(tempb);
                }
            }
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_mul(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb;
    MathsNode **tempb;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&bb,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = mulGenVar(aa,bb);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }

        else
        {
            if ( isReal(aa) )
            {
                /* the first arg cannot be zero */

                tempb = _this->SubFunctions;

                (tempb[1])->multiplier = mulGenVar(_this->multiplier,(tempb[1])->multiplier);
                (tempb[1])->multiplier = mulGenVar(aa,(tempb[1])->multiplier);

                *_this = *(tempb[1]);

                delEqn(tempb[0]);
                free(tempb[1]);

                free(tempb);
            }

            else
            {
                if ( isReal(bb) )
                {
                    /* the first arg cannot be zero */

                    tempb = _this->SubFunctions;

                    (tempb[0])->multiplier = mulGenVar(_this->multiplier,(tempb[0])->multiplier);
                    (tempb[0])->multiplier = mulGenVar(bb,(tempb[0])->multiplier);

                    *_this = *(tempb[0]);

                    free(tempb[0]);
                    delEqn(tempb[1]);

                    free(tempb);
                }
            }
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_div(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb;
    MathsNode **temp;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&bb,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = divGenVar(aa,bb);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }

        else
        {
            if ( isReal(bb) )
            {
                temp = _this->SubFunctions;

                (temp[0])->multiplier = mulGenVar(_this->multiplier,(temp[0])->multiplier);
                (temp[0])->multiplier = divGenVar((temp[0])->multiplier,bb);

                *_this = *(temp[0]);

                free(temp[0]);
                delEqn(temp[1]);

                free(temp);
            }
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_idiv(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&bb,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = idivGenVar(aa,bb);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_fmod(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&bb,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = modGenVar(aa,bb);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_eval(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar tempb;
    long m;
    GenVar mm,ii;
    MathsNode **xemp;
    GenVar dummy;
    int q;

    if ( ( q = locevaluateEqn(&mm,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( isFiniteInteger(mm) )
    {
        if ( simplify_onthefly )
        {
            /*
              OK, m is just a simple integer, so we can remove the eval
              function in this case.
            */

            xemp = _this->SubFunctions;

            (xemp[2])->multiplier = mulGenVar(_this->multiplier,(xemp[2])->multiplier);

            *_this = *(xemp[2]);

            locsubVar(_this,xemp[1],0,getInteger(mm),0);

            delEqn(xemp[0]);
            delEqn(xemp[1]);
            free(xemp[2]);

            free(xemp);

            /*
              we have replaced var(0,m) with i in x and replaced the old
              expression with x (after substitution).  This may have
              resulted in further possible simplification, so we transfer
              control to the relevant function to complete any further
              simplifications.
            */

            return ((FunctionTable[_this->CombinerFnNumber]).CalcOperation)(result,_this,GlobInput,touch_all,report_change,simplify_onthefly);
        }

        else
        {
            /*
              We need i, so let's get it.
            */

            if ( ( q = locevaluateEqn(&ii,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            m = getInteger(mm);

            tempb = get_loc_var(m,GlobInput);

            if ( touch_loc_var(m,GlobInput) == 0 ) /* this will return -1 oob */
            {
                set_loc_var(m,ii,GlobInput);

                if ( ( q = locevaluateEqn(result,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return q;
                }
            }

            else
            {
                (*result) = tempb;

                if ( touch_all )
                {
                    if ( ( q = locevaluateEqn(&dummy,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                    {
                        return q;
                    }
                }
            }

            restore_loc_var(m,tempb,GlobInput);
        }
    }

    else
    {
        if ( mm.DataType == DATA_IS_UNCALCED )
        {
            /*
               Big caveat here: The default (ie. just after initialisation)
               value for all local variables is DATA_IS_UNCALCED.  Thus, if
               the local variable is used here then it will appear as such,
               meaning that the result here will be uncalced also, which is
               what we want.  If the local variable is unnecessary then the
               result here will be a number (or error) which is independent
               of the local variable.  Thus this is a valid statement,
               *ASSUMING THAT THERE IS NO RECURSIVE REUSE HAPPENING*.
            */

            if ( ( q = locevaluateEqn(&ii,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }

            if ( ( q = locevaluateEqn(result,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return q;
            }
        }

        else
        {
            return EERR_INVALIDVAR;
        }
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_sum(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar temp;
    GenVar tempb;
    long m,i,j,k,q;
    double di,dj,dk,dq;
    GenVar mm,ii,jj,kk,zzz;
    MathsNode **xemp;
    MathsNode *xempb;
    GenVar dummy;
    int s;

    if ( ( s = locevaluateEqn(&mm,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    if ( ( s = locevaluateEqn(&ii,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    if ( ( s = locevaluateEqn(&jj,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    if ( ( s = locevaluateEqn(&kk,(_this->SubFunctions)[3],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    zzz = ResultUncalced();

    if ( isFiniteReal(ii) && isFiniteReal(jj) && simplify_onthefly )
    {
        if ( isaequalb(ii,jj) )
        {
            /*
              In this case, mm = ii = jj, so we can replace with eval and
              let that do the simplification.
            */

            xemp = _this->SubFunctions;

            xempb = copyEquation(xsumprod_simplifier);

            xempb->multiplier = mulGenVar(xempb->multiplier,_this->multiplier);

            *_this = *xempb;

            free(xempb);

            locsubVar(_this,xemp[0],-1,-1,0);
            locsubVar(_this,xemp[1],-2,-2,0);
            locsubVar(_this,xemp[4],-3,-3,0);

            delEqn(xemp[0]);
            delEqn(xemp[1]);
            delEqn(xemp[2]);
            delEqn(xemp[3]);
            delEqn(xemp[4]);

            free(xemp);

            /*
              Transfer control for further possible simplification.
            */

            return ((FunctionTable[_this->CombinerFnNumber]).CalcOperation)(result,_this,GlobInput,touch_all,report_change,simplify_onthefly);
        }
    }

    if ( isFiniteInteger(mm) && isFiniteReal(ii) && isFiniteReal(jj) && isFiniteReal(kk) )
    {
        if ( simplify_onthefly && 1 ) /* FIXME: 0 will always evaluate false,
                                         thus temporarily excising this
                                         code.  But for reasons unknown this
                                         will break differentiation, so
                                         temporarily re-allow to make this
                                         work (for now). */
        {
            /*
               Special case: the entire sum is well defined, so we just
               unroll the summation recursively to remove the summation.

               Note: we cannot have ii == jj here, as that case have already
               been dealt with.

               FIXME: This may fail with reals, as the end of the summation
                      may not correctly line up (due to numerical errors),
                      making the following test trip incorrectly.
            */

            if ( (  isalessb(jj,ii) && getReal(kk) >= 0 ) ||
                 ( !isalessb(jj,ii) && getReal(kk) <= 0 )    )
            {
                (*result) = ResultError(GVAR_EILLINTEG);
            }

            else
            {
                xemp = _this->SubFunctions;

                xempb = copyEquation(xsum_recurser);

                xempb->multiplier = mulGenVar(xempb->multiplier,_this->multiplier);

                *_this = *xempb;

                free(xempb);

                locsubVar(_this,xemp[0],-1,-1,0);
                locsubVar(_this,xemp[1],-2,-2,0);
                locsubVar(_this,xemp[2],-3,-3,0);
                locsubVar(_this,xemp[3],-4,-4,0);
                locsubVar(_this,xemp[4],-5,-5,0);

                delEqn(xemp[0]);
                delEqn(xemp[1]);
                delEqn(xemp[2]);
                delEqn(xemp[3]);
                delEqn(xemp[4]);

                free(xemp);

                /*
                  Transfer control for further possible simplification.
                */

                return ((FunctionTable[_this->CombinerFnNumber]).CalcOperation)(result,_this,GlobInput,touch_all,report_change,simplify_onthefly);
            }
        }

        else
        {
            m = getInteger(mm);

            tempb = get_loc_var(m,GlobInput);

            if ( isFiniteInteger(ii) && isFiniteInteger(jj) && isFiniteInteger(kk) )
            {
                i = getInteger(ii);
                j = getInteger(jj);
                k = getInteger(kk);

                if ( j < i )
                {
                    q = j;
                    j = i;
                    i = q;
                    k = -k;
                }

                if ( k <= 0 )
                {
                    (*result) = ResultError(GVAR_EILLINTEG);

                    if ( touch_all )
                    {
                        if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                        {
                            return s;
                        }
                    }
                }

                else
                {
                    if ( touch_loc_var(m,GlobInput) == 0 ) /* this will return -1 oob */
                    {
                        (*result) = ResultZero();

                        for ( q = i ; q <= j ; q = q + k )
                        {
                            temp = ResultInteger(q);
                            set_loc_var(m,temp,GlobInput);
                            if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                            {
                                return s;
                            }

                            (*result) = addGenVar((*result),dummy);

                            if ( isError((*result)) )
                            {
                                break;
                            }
                        }
                    }

                    else
                    {
                        (*result) = tempb;

                        if ( touch_all )
                        {
                            if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                            {
                                return s;
                            }
                        }
                    }
                }
            }

            else
            {
                di = getReal(ii);
                dj = getReal(jj);
                dk = getReal(kk);

                if ( dj < di )
                {
                    dq = dj;
                    dj = di;
                    di = dq;
                    dk = -dk;
                }

                if ( dk <= 0 )
                {
                    (*result) = ResultError(GVAR_EILLINTEG);

                    if ( touch_all )
                    {
                        if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                        {
                            return s;
                        }
                    }
                }

                else
                {
                    if ( touch_loc_var(m,GlobInput) == 0 ) /* this will return -1 oob */
                    {
                        (*result) = ResultZero();

                        for ( dq = di ; dq <= dj ; dq = dq + dk )
                        {
                            temp = ResultReal(dq);
                            set_loc_var(m,temp,GlobInput);
                            if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                            {
                                return s;
                            }

                            (*result) = addGenVar((*result),dummy);

                            if ( isError((*result)) )
                            {
                                break;
                            }
                        }
                    }

                    else
                    {
                        (*result) = tempb;

                        if ( touch_all )
                        {
                            if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                            {
                                return s;
                            }
                        }
                    }
                }
            }

            restore_loc_var(m,tempb,GlobInput);
        }
    }

    else
    {
        if ( ( ( mm.DataType == DATA_IS_UNCALCED ) || isFiniteInteger(mm) ) &&
             ( ( ii.DataType == DATA_IS_UNCALCED ) || isReal(ii)          ) &&
             ( ( jj.DataType == DATA_IS_UNCALCED ) || isReal(jj)          ) &&
             ( ( kk.DataType == DATA_IS_UNCALCED ) || isReal(kk)          )    )
        {
            /*
               Big caveat here: The default (ie. just after initialisation)
               value for all local variables is DATA_IS_UNCALCED.  Thus, if
               the local variable is used here then it will appear as such,
               meaning that the result here will be uncalced also, which is
               what we want.  If the local variable is unnecessary then the
               result here will be a number (or error) which is independent
               of the local variable.  Thus this is a valid statement,
               *ASSUMING THAT THERE IS NO RECURSIVE REUSE HAPPENING*.
            */

            if ( ( s = locevaluateEqn(&zzz,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return s;
            }

            (*result) = ResultUncalced();
        }

        else
        {
            return EERR_INVALIDVAR;
        }
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_prod(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar temp;
    GenVar tempb;
    long m,i,j,k,q;
    double di,dj,dk,dq;
    GenVar mm,ii,jj,kk,zzz;
    MathsNode **xemp;
    MathsNode *xempb;
    GenVar dummy;
    int s;

    if ( ( s = locevaluateEqn(&mm,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    if ( ( s = locevaluateEqn(&ii,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    if ( ( s = locevaluateEqn(&jj,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    if ( ( s = locevaluateEqn(&kk,(_this->SubFunctions)[3],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    zzz = ResultUncalced();

    if ( isFiniteReal(ii) && isFiniteReal(jj) && simplify_onthefly )
    {
        if ( isaequalb(ii,jj) )
        {
            /*
              In this case, mm = ii = jj, so we can replace with eval and
              let that do the simplification.
            */

            xemp = _this->SubFunctions;

            xempb = copyEquation(xsumprod_simplifier);

            xempb->multiplier = mulGenVar(xempb->multiplier,_this->multiplier);

            *_this = *xempb;

            free(xempb);

            locsubVar(_this,xemp[0],-1,-1,0);
            locsubVar(_this,xemp[1],-2,-2,0);
            locsubVar(_this,xemp[4],-3,-3,0);

            delEqn(xemp[0]);
            delEqn(xemp[1]);
            delEqn(xemp[2]);
            delEqn(xemp[3]);
            delEqn(xemp[4]);

            free(xemp);

            /*
              Transfer control for further possible simplification.
            */

            return ((FunctionTable[_this->CombinerFnNumber]).CalcOperation)(result,_this,GlobInput,touch_all,report_change,simplify_onthefly);
        }
    }

    if ( isFiniteInteger(mm) && isFiniteReal(ii) && isFiniteReal(jj) && isFiniteReal(kk) )
    {
        if ( simplify_onthefly && 1 ) /* FIXME: 0 will always evaluate false,
                                         thus temporarily excising this
                                         code.  But for reasons unknown this
                                         will break differentiation, so
                                         temporarily re-allow to make this
                                         work (for now). */
        {
            /*
               Special case: the entire product is well defined, so we just
               unroll it recursively to remove the production.

               Note: we cannot have ii == jj here, as that case have already
               been dealt with.

               FIXME: This may fail with reals, as the end of the prodation
                      may not correctly line up (due to numerical errors),
                      making the following test trip incorrectly.
            */

            if ( (  isalessb(jj,ii) && getReal(kk) >= 0 ) ||
                 ( !isalessb(jj,ii) && getReal(kk) <= 0 )    )
            {
                (*result) = ResultError(GVAR_EILLINTEG);
            }

            else
            {
                xemp = _this->SubFunctions;

                xempb = copyEquation(xprod_recurser);

                xempb->multiplier = mulGenVar(xempb->multiplier,_this->multiplier);

                *_this = *xempb;

                free(xempb);

                locsubVar(_this,xemp[0],-1,-1,0);
                locsubVar(_this,xemp[1],-2,-2,0);
                locsubVar(_this,xemp[2],-3,-3,0);
                locsubVar(_this,xemp[3],-4,-4,0);
                locsubVar(_this,xemp[4],-5,-5,0);

                delEqn(xemp[0]);
                delEqn(xemp[1]);
                delEqn(xemp[2]);
                delEqn(xemp[3]);
                delEqn(xemp[4]);

                free(xemp);

                /*
                  Transfer control for further possible simplification.
                */

                return ((FunctionTable[_this->CombinerFnNumber]).CalcOperation)(result,_this,GlobInput,touch_all,report_change,simplify_onthefly);
            }
        }

        else
        {
            m = getInteger(mm);

            tempb = get_loc_var(m,GlobInput);

            if ( isFiniteInteger(ii) && isFiniteInteger(jj) && isFiniteInteger(kk) )
            {
                i = getInteger(ii);
                j = getInteger(jj);
                k = getInteger(kk);

                if ( j < i )
                {
                    q = j;
                    j = i;
                    i = q;
                    k = -k;
                }

                if ( k <= 0 )
                {
                    (*result) = ResultError(GVAR_EILLINTEG);

                    if ( touch_all )
                    {
                        if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                        {
                            return s;
                        }
                    }
                }

                else
                {
                    if ( touch_loc_var(m,GlobInput) == 0 )
                    {
                        (*result) = ResultOne();

                        for ( q = i ; q <= j ; q = q + k )
                        {
                            temp = ResultInteger(q);
                            set_loc_var(m,temp,GlobInput);
                            if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                            {
                                return s;
                            }

                            (*result) = mulGenVar((*result),dummy);

                            if ( isError((*result)) )
                            {
                                break;
                            }
                        }
                    }

                    else
                    {
                        (*result) = tempb;

                        if ( touch_all )
                        {
                            if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                            {
                                return s;
                            }
                        }
                    }
                }
            }

            else
            {
                di = getInteger(ii);
                dj = getInteger(jj);
                dk = getInteger(kk);

                if ( dj < di )
                {
                    dq = dj;
                    dj = di;
                    di = dq;
                    dk = -dk;
                }

                if ( dk <= 0 )
                {
                    (*result) = ResultError(GVAR_EILLINTEG);

                    if ( touch_all )
                    {
                        if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                        {
                            return s;
                        }
                    }
                }

                else
                {
                    if ( touch_loc_var(m,GlobInput) == 0 )
                    {
                        (*result) = ResultOne();

                        for ( dq = di ; dq <= dj ; dq = dq + dk )
                        {
                            temp = ResultReal(dq);
                            set_loc_var(m,temp,GlobInput);
                            if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                            {
                                return s;
                            }

                            (*result) = mulGenVar((*result),dummy);

                            if ( isError((*result)) )
                            {
                                break;
                            }
                        }
                    }

                    else
                    {
                        (*result) = tempb;

                        if ( touch_all )
                        {
                            if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                            {
                                return s;
                            }
                        }
                    }
                }
            }

            restore_loc_var(m,tempb,GlobInput);
        }
    }

    else
    {
        if ( ( ( mm.DataType == DATA_IS_UNCALCED ) || isFiniteInteger(mm) ) &&
             ( ( ii.DataType == DATA_IS_UNCALCED ) || isReal(ii)          ) &&
             ( ( jj.DataType == DATA_IS_UNCALCED ) || isReal(jj)          ) &&
             ( ( kk.DataType == DATA_IS_UNCALCED ) || isReal(kk)          )    )
        {
            /*
               Big caveat here: The default (ie. just after initialisation)
               value for all local variables is DATA_IS_UNCALCED.  Thus, if
               the local variable is used here then it will appear as such,
               meaning that the result here will be uncalced also, which is
               what we want.  If the local variable is unnecessary then the
               result here will be a number (or error) which is independent
               of the local variable.  Thus this is a valid statement,
               *ASSUMING THAT THERE IS NO RECURSIVE REUSE HAPPENING*.
            */

            if ( ( s = locevaluateEqn(&zzz,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
            {
                return s;
            }

            (*result) = ResultUncalced();
        }

        else
        {
            return EERR_INVALIDVAR;
        }
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_integ(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar temp;
    GenVar tempb;
    long m;
    double i,j,k,q,r;
    GenVar mm,ii,jj,kk,zzz;
    MathsNode **xemp;
    MathsNode *xempb;
    GenVar dummy;
    int s;

    if ( ( s = locevaluateEqn(&mm,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    if ( ( s = locevaluateEqn(&ii,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    if ( ( s = locevaluateEqn(&jj,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    if ( ( s = locevaluateEqn(&kk,(_this->SubFunctions)[3],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    if ( simplify_onthefly )
    {
        if ( isFiniteReal(ii) && isFiniteReal(jj) )
        {
            if ( isaequalb(ii,jj) )
            {
                /*
                  The result is 0 (note flotilla of assumptions here).
                */

                xemp = _this->SubFunctions;

                *_this = *xglobal_zero;

                delEqn(xemp[0]);
                delEqn(xemp[1]);
                delEqn(xemp[2]);
                delEqn(xemp[3]);
                delEqn(xemp[4]);

                free(xemp);

                (*result) = ResultZero();

                return 0;
            }
        }

        if ( ( s = locevaluateEqn(&zzz,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
        {
            return s;
        }

        if ( isReal(zzz) )
        {
            xemp = _this->SubFunctions;

            xempb = copyEquation(xinteg_simplifier);

            xempb->multiplier = mulGenVar(xempb->multiplier,_this->multiplier);

            *_this = *xempb;

            free(xempb);

            locsubVar(_this,xemp[4],-1,-1,0);
            locsubVar(_this,xemp[2],-2,-2,0);
            locsubVar(_this,xemp[1],-3,-3,0);

            delEqn(xemp[0]);
            delEqn(xemp[1]);
            delEqn(xemp[2]);
            delEqn(xemp[3]);
            delEqn(xemp[4]);

            free(xemp);

            /*
              Transfer control for further possible simplification.
            */

            return ((FunctionTable[_this->CombinerFnNumber]).CalcOperation)(result,_this,GlobInput,touch_all,report_change,simplify_onthefly);
        }
    }

    else
    {
        zzz = ResultUncalced();
    }

    if ( isFiniteInteger(mm) && isFiniteReal(ii) && isFiniteReal(jj) && isFiniteReal(kk) )
    {
        m = getInteger(mm);
        i = getReal(ii);
        j = getReal(jj);
        k = getReal(kk);

        tempb = get_loc_var(m,GlobInput);

        if ( j < i )
        {
            q = j;
            j = i;
            i = q;

            if ( k == 0 )
            {
                k = DEFAULT_INTEGRAL_STEP;
            }

            else
            {
                k = -k;
            }
        }

        else
        {
            if ( k == 0 )
            {
                k = DEFAULT_INTEGRAL_STEP;
            }
        }

        if ( i == j )
        {
            (*result) = ResultZero();

            if ( touch_all && !simplify_onthefly )
            {
                if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return s;
                }
            }
        }

        else
        {
            if ( k <= 0 )
            {
                (*result) = ResultError(GVAR_EILLINTEG);

                if ( touch_all && !simplify_onthefly )
                {
                    if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                    {
                        return s;
                    }
                }
            }

            else
            {
                kk = ResultReal(k*(j-i));

                if ( touch_loc_var(m,GlobInput) == 0 )
                {
                    (*result) = ResultZero();

                    for ( q = 0 ; q <= 1 ; q = q + k )
                    {
                        r = q;
                        r *= (j-i);
                        r += i;

                        temp = ResultReal(r);
                        set_loc_var(m,temp,GlobInput);
                        if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                        {
                            return s;
                        }

                        temp = mulGenVar(kk,dummy);
                        (*result) = addGenVar((*result),temp);

                        if ( isError((*result)) )
                        {
                            break;
                        }
                    }
                }

                else
                {
                    (*result) = tempb;

                    if ( touch_all && !simplify_onthefly )
                    {
                        if ( ( s = locevaluateEqn(&dummy,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                        {
                            return s;
                        }
                    }
                }
            }
        }

        restore_loc_var(m,tempb,GlobInput);
    }

    else
    {
        if ( ( ( mm.DataType == DATA_IS_UNCALCED ) || isFiniteInteger(mm) ) &&
             ( ( ii.DataType == DATA_IS_UNCALCED ) || isReal(ii)          ) &&
             ( ( jj.DataType == DATA_IS_UNCALCED ) || isReal(jj)          ) &&
             ( ( kk.DataType == DATA_IS_UNCALCED ) || isReal(kk)          )    )
        {
            /*
               Big caveat here: The default (ie. just after initialisation)
               value for all local variables is DATA_IS_UNCALCED.  Thus, if
               the local variable is used here then it will appear as such,
               meaning that the result here will be uncalced also, which is
               what we want.  If the local variable is unnecessary then the
               result here will be a number (or error) which is independent
               of the local variable.  Thus this is a valid statement,
               *ASSUMING THAT THERE IS NO RECURSIVE REUSE HAPPENING*.
            */

            if ( !simplify_onthefly && touch_all )
            {
                /*
                  no need to do this twice.
                */

                if ( ( s = locevaluateEqn(&zzz,(_this->SubFunctions)[4],GlobInput,touch_all,report_change,simplify_onthefly) ) )
                {
                    return s;
                }
            }

            (*result) = ResultUncalced();
        }

        else
        {
            return EERR_INVALIDVAR;
        }
    }

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_max(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&bb,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_max(aa,bb);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_min(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&bb,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_min(aa,bb);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_abs(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_abs(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_sgn(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_sgn(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_rint(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_rint(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_ceil(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_ceil(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_floor(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_floor(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_inv(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = invGenVar(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_neg(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = negGenVar(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_sin(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_sin(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_cos(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_cos(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_tan(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_tan(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_cosec(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_cosec(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_sec(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_sec(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_cot(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_cot(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_asin(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_asin(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_acos(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_acos(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_atan(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_atan(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_acosec(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_acosec(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_asec(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_asec(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_acot(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_acot(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_sinh(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_sinh(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_cosh(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_cosh(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_tanh(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_tanh(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_cosech(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_cosech(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_sech(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_sech(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_coth(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_coth(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_asinh(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_asinh(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_acosh(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_acosh(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_atanh(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_atanh(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_acosech(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_acosech(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_asech(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_asech(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_acoth(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_acoth(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_sinc(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_sinc(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_gd(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_gd(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_agd(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_agd(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_exp(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_exp(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_tenup(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_tenup(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_log(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_log(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_logX(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_logX(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_pow(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar xx,yy;
    MathsNode **temp;
    int q;

    if ( ( q = locevaluateEqn(&xx,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&yy,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_pow(xx,yy);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }

        else
        {
            if ( yy.DataType == DATA_IS_ONE )
            {
                temp = _this->SubFunctions;

                _this->CombinerFnNumber = (temp[0])->CombinerFnNumber;
                _this->multiplier       = mulGenVar((temp[0])->multiplier,_this->multiplier);
                _this->SubFunctions     = (temp[0])->SubFunctions;

                free(temp[0]);
                delEqn(temp[1]);

                free(temp);
            }
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_sqrt(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_sqrt(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_gamma(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_gamma(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_lngamma(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_lngamma(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_normDistr(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;
    
    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_normDistr(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_polyDistr(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&bb,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_polyDistr(aa,bb);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_erf(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_erf(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_erfc(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_erfc(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_perm(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&bb,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_perm(aa,bb);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_comb(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&bb,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_comb(aa,bb);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_fact(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_fact(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_randn(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    double xx,yy;
    GenVar x,y;
    int q;

    RAND_INIT_SEQUENCE;

    if ( report_change )
    {
        (*result) = ResultUncalced();

        return 0;
    }

    if ( ( q = locevaluateEqn(&x,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&y,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    switch ( x.DataType )
    {
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        case DATA_IS_ONE:
        case DATA_IS_ZERO:
        {
            break;
        }

        case DATA_IS_POS_INFTY:
        {
            x = ResultReal(DBL_MAX);

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            x = ResultReal(-DBL_MAX);

            break;
        }

        case DATA_IS_UNCALCED:
        {
            (*result) = x;

            return 0;

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_ERROR:
        default:
        {
            (*result) = addGenVar(x,y);

            return 0;

            break;
        }
    }

    switch ( y.DataType )
    {
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        case DATA_IS_ONE:
        case DATA_IS_ZERO:
        {
            break;
        }

        case DATA_IS_POS_INFTY:
        {
            y = ResultReal(DBL_MAX);

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            y = ResultReal(-DBL_MAX);

            break;
        }

        case DATA_IS_UNCALCED:
        {
            (*result) = y;

            return 0;

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_ERROR:
        default:
        {
            (*result) = addGenVar(x,y);

            return 0;

            break;
        }
    }

    xx = getReal(x);
    yy = getReal(y);

    (*result) = ResultReal(RandomDouble(xx,yy));

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_irand(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    long xx,yy;
    GenVar x,y;
    int q;

    RAND_INIT_SEQUENCE;

    if ( report_change )
    {
        (*result) = ResultUncalced();

        return 0;
    }

    if ( ( q = locevaluateEqn(&x,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&y,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    switch ( x.DataType )
    {
        case DATA_IS_INTEGER:
        {
            break;
        }

        case DATA_IS_REAL:
        {
            x = Gen_ceil(x);

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_ZERO:
        {
            break;
        }

        case DATA_IS_POS_INFTY:
        {
            x = ResultInteger(LONG_MAX);

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            x = ResultInteger(LONG_MIN);

            break;
        }

        case DATA_IS_UNCALCED:
        {
            (*result) = x;

            return 0;

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_ERROR:
        default:
        {
            (*result) = addGenVar(x,y);

            return 0;

            break;
        }
    }

    switch ( y.DataType )
    {
        case DATA_IS_INTEGER:
        {
            break;
        }

        case DATA_IS_REAL:
        {
            y = Gen_floor(y);

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_ZERO:
        {
            break;
        }

        case DATA_IS_POS_INFTY:
        {
            y = ResultInteger(LONG_MAX);

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            y = ResultInteger(LONG_MIN);

            break;
        }

        case DATA_IS_UNCALCED:
        {
            (*result) = y;

            return 0;

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_ERROR:
        default:
        {
            (*result) = addGenVar(x,y);

            return 0;

            break;
        }
    }

    xx = getInteger(x);
    yy = getInteger(y);

    (*result) = ResultReal(RandomInt(xx,yy));

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_grand(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar mean,stddev;
    double q;
    int s;

    RAND_INIT_SEQUENCE;

    if ( report_change )
    {
        (*result) = ResultUncalced();

        return 0;
    }

    if ( ( s = locevaluateEqn(&mean,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    if ( ( s = locevaluateEqn(&stddev,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    switch ( stddev.DataType )
    {
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        case DATA_IS_ONE:
        {
            q = getReal(stddev);

            if ( q > 0 )
            {
                (*result) = ResultReal(RandomGaussian(0,q));
            }

            else
            {
                (*result) = ResultError(GVAR_GSL_EDOM);
            }

            break;
        }

        case DATA_IS_ZERO:
        {
            (*result) = ResultZero();

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            (*result) = ResultReal(RandomDouble(-DBL_MAX,DBL_MAX));

            break;
        }

        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            (*result) = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_IS_UNCALCED:
        {
            (*result) = stddev;

            return 0;

            break;
        }

        case DATA_IS_ERROR:
        default:
        {
            (*result) = stddev;

            break;
        }
    }

    (*result) = addGenVar((*result),mean);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_prand(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar mean,stddev;
    double q;
    GenVar order;
    long p;
    int s;

    RAND_INIT_SEQUENCE;

    if ( report_change )
    {
        (*result) = ResultUncalced();

        return 0;
    }

    if ( ( s = locevaluateEqn(&mean,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    if ( ( s = locevaluateEqn(&stddev,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    if ( ( s = locevaluateEqn(&order,(_this->SubFunctions)[2],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return s;
    }

    if ( isFiniteInteger(order) )
    {
        switch ( stddev.DataType )
        {
            case DATA_IS_INTEGER:
            case DATA_IS_REAL:
            case DATA_IS_ONE:
            {
                q = getReal(stddev);
                p = getInteger(order);

		if ( ( q > 0 ) && ( isFiniteInteger(order) ) && ( p > 0 ) && ( p < 12 ) ) // FIXME: need to be able to do more than 11th order poly noise
                {
                    (*result) = ResultReal(RandomPolynomial(0,q,p));
                }

                else
                {
                    (*result) = ResultError(GVAR_GSL_EDOM);
                }

                break;
            }

            case DATA_IS_ZERO:
            {
                (*result) = ResultZero();

                break;
            }

            case DATA_IS_POS_INFTY:
            {
                (*result) = ResultReal(RandomDouble(-DBL_MAX,DBL_MAX));

                break;
            }

            case DATA_IS_NEG_INFTY:
            case DATA_IS_FIN_INDET:
            case DATA_IS_INF_INDET:
            {
                (*result) = ResultError(GVAR_GSL_EDOM);

                break;
            }

            case DATA_IS_UNCALCED:
            {
                (*result) = stddev;

                return 0;

                break;
            }

            case DATA_IS_ERROR:
            default:
            {
                (*result) = stddev;

                break;
            }
        }
    }

    else
    {
        switch ( order.DataType )
        {
            case DATA_IS_REAL:
            case DATA_IS_INTEGER:
            case DATA_IS_ONE:
            case DATA_IS_ZERO:
            case DATA_IS_POS_INFTY:
            case DATA_IS_NEG_INFTY:
            case DATA_IS_FIN_INDET:
            case DATA_IS_INF_INDET:
            {
                (*result) = ResultError(GVAR_GSL_EDOM);

                break;
            }

            case DATA_IS_UNCALCED:
            {
                (*result) = order;

                return 0;

                break;
            }

            case DATA_IS_ERROR:
            default:
            {
                (*result) = order;

                break;
            }
        }
    }

    (*result) = addGenVar((*result),mean);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_gami(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&bb,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_gami(aa,bb);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_gamic(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&bb,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_gamic(aa,bb);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_psi(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_psi(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_psi_n(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa,bb;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    if ( ( q = locevaluateEqn(&bb,(_this->SubFunctions)[1],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_psi_n(aa,bb);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}

int Calc_dawson(GenVar *result, MathsNode *_this, GenArg GlobInput, int touch_all, int report_change, int simplify_onthefly)
{
    GenVar aa;
    int q;

    if ( ( q = locevaluateEqn(&aa,(_this->SubFunctions)[0],GlobInput,touch_all,report_change,simplify_onthefly) ) )
    {
        return q;
    }

    (*result) = Gen_dawson(aa);

    if ( simplify_onthefly )
    {
        if ( (*result).DataType != DATA_IS_UNCALCED )
        {
            make_number(_this,result);
        }
    }

    return 0;

    REMOVE_GCC_WARNING_CALC;
}





/*
   Stuff borrowed from elsewhere
   NOTE: I did not write this code.
*/

/*
   This Random Number Generator is based on the algorithm in a FORTRAN
   version published by George Marsaglia and Arif Zaman, Florida State
   University; ref.: see original comments below.
   At the fhw (Fachhochschule Wiesbaden, W.Germany), Dept. of Computer
   Science, we have written sources in further languages (C, Modula-2
   Turbo-Pascal(3.0, 5.0), Basic and Ada) to get exactly the same test
   results compared with the original FORTRAN version.
   April 1989
   Karl-L. Noell <NOELL@DWIFH1.BITNET>
      and  Helmut  Weber <WEBER@DWIFH1.BITNET>

   This random number generator originally appeared in "Toward a Universal
   Random Number Generator" by George Marsaglia and Arif Zaman.
   Florida State University Report: FSU-SCRI-87-50 (1987)
   It was later modified by F. James and published in "A Review of Pseudo-
   random Number Generators"
   THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
   (However, a newly discovered technique can yield
   a period of 10^600. But that is still in the development stage.)
   It passes ALL of the tests for random number generators and has a period
   of 2^144, is completely portable (gives bit identical results on all
   machines with at least 24-bit mantissas in the floating point
   representation).
   The algorithm is a combination of a Fibonacci sequence (with lags of 97
   and 33, and operation "subtraction plus one, modulo one") and an
   "arithmetic sequence" (using subtraction).

   Use IJ = 1802 & KL = 9373 to test the random number generator. The
   subroutine RANMAR should be used to generate 20000 random numbers.
   Then display the next six random numbers generated multiplied by 4096*4096
   If the random number generator is working properly, the random numbers
   should be:
           6533892.0  14220222.0  7275067.0
           6172232.0  8354498.0   10633180.0
*/

/* Globals */
#define FALSE   0
#define TRUE    1

double u[97],c,cd,cm;
int i97,j97;
int test = FALSE;

/*
   This is the initialization routine for the random number generator.
   NOTE: The seed variables can have values between:    0 <= IJ <= 31328
                                                        0 <= KL <= 30081
   The random number sequences created by these two seeds are of sufficient
   length to complete an entire calculation with. For example, if sveral
   different groups are working on different parts of the same calculation,
   each group could be assigned its own IJ seed. This would leave each group
   with 30000 choices for the second seed. That is to say, this random
   number generator can create 900 million different subsequences -- with
   each subsequence having a length of approximately 10^30.
*/
void RandomInitialise(int ij,int kl)
{
   double s,t;
   int ii,i,j,k,l,jj,m;

   /*
      Handle the seed range errors
         First random number seed must be between 0 and 31328
         Second seed must have a value between 0 and 30081
   */
   if (ij < 0 || ij > 31328 || kl < 0 || kl > 30081) {
      ij = 1802;
      kl = 9373;
   }

   i = (ij / 177) % 177 + 2;
   j = (ij % 177)       + 2;
   k = (kl / 169) % 178 + 1;
   l = (kl % 169);

   for (ii=0; ii<97; ii++) {
      s = 0.0;
      t = 0.5;
      for (jj=0; jj<24; jj++) {
         m = (((i * j) % 179) * k) % 179;
         i = j;
         j = k;
         k = m;
         l = (53 * l + 1) % 169;
         if (((l * m % 64)) >= 32)
            s += t;
         t *= 0.5;
      }
      u[ii] = s;
   }

   c    = 362436.0 / 16777216.0;
   cd   = 7654321.0 / 16777216.0;
   cm   = 16777213.0 / 16777216.0;
   i97  = 97;
   j97  = 33;
   test = TRUE;
}

/* 
   This is the random number generator proposed by George Marsaglia in
   Florida State University Report: FSU-SCRI-87-50
*/
double RandomUniform(void)
{
   double uni;

   /* Make sure the initialisation routine has been called */
   if (!test) 
      RandomInitialise(1802,9373);

   uni = u[i97-1] - u[j97-1];
   if (uni <= 0.0)
      uni++;
   u[i97-1] = uni;
   i97--;
   if (i97 == 0)
      i97 = 97;
   j97--;
   if (j97 == 0)
      j97 = 97;
   c -= cd;
   if (c < 0.0)
      c += cm;
   uni -= c;
   if (uni < 0.0)
      uni++;

   return(uni);
}

/*
  ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
  THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
  VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
  The function returns a normally distributed pseudo-random number
  with a given mean and standard devaiation.  Calls are made to a
  function subprogram which must return independent random
  numbers uniform in the interval (0,1).
  The algorithm uses the ratio of uniforms method of A.J. Kinderman
  and J.F. Monahan augmented with quadratic bounding curves.
*/
double RandomGaussian(double mean,double stddev)
{
    /*
       FIXME: this all looks very nice and fancy, but it DOES NOT
       GENERATE NUMBERS WITH APPROPRIATE MEAN AND VARIANCE.
       It just doesn't.  So I've replaced the whole schebang with
       a simple call to the polynomial noise one for degree 2.
    */

   return RandomPolynomial(mean,stddev,2);

   double  q,u,v,x,y;

   /*  
      Generate P = (u,v) uniform in rect. enclosing acceptance region 
      Make sure that any random numbers <= 0 are rejected, since
      gaussian() requires uniforms > 0, but RandomUniform() delivers >= 0.
   */
   do {
      u = RandomUniform();
      v = RandomUniform();
      if (u <= 0.0 || v <= 0.0) {
          u = 1.0;
          v = 1.0;
      }
      v = 1.7156 * (v - 0.5);

      /*  Evaluate the quadratic form */
      x = u - 0.449871;
      y = abs(v) + 0.386595;
      q = x * x + y * (0.19600 * y - 0.25472 * x);

      /* Accept P if inside inner ellipse */
      if (q < 0.27597)
         break;

      /*  Reject P if outside outer ellipse, or outside acceptance region */
    } while ((q > 0.27846) || (v * v > -4.0 * log(u) * u * u));

    /*  Return ratio of P's coordinates as the normal deviate */
    return (mean + stddev * v / u);
}





/*
My adaption of RandomGaussian to give numbers from a polynomial
distribution.  Actually, it only approximates the polynomial distribution
to the 95% point.  Infinite loop is possible, but (I hope) unlikely.
*/

double RandomPolynomial(double mean, double stddev, int order)
{
    double u,v,x,y,z,cp,cp_prime,edge,forder;
    double tempa;
    double tempb;

    /*
       This code works by generating a random number within a certain
       range and then randomly choosing whether to accept this number
       based on the probability that a number from the actual distribution
       would fall at this point.

       The uniform distribution range is selected to cover the 99% percentile,
       but for practical reasons must be expanded for regions that are
       too narrow (lower order polynomial distributions) and contracted
       for regions that are too wide (higher order polynomial
       distributions).  The point of expansion is largely arbitrary (and
       has no real impact).  The point of contraction comes from practical
       observation of when the variance of the result diverges (11th order
       if no contraction is present).

       NB: This will not be terribly accurate for 12th or higher order.
    */

    forder = (double) order;

    tempa = 1/forder;
    tempb = 3/forder;

    x = gsl_sf_gamma(tempa);
    y = gsl_sf_gamma(tempb);
    z = sqrt(y/x);

    cp       = 0.5*forder*z/x;
    cp_prime = pow(z,forder);

    edge = log(cp/0.05)/cp_prime;

    if ( edge <= 100 )
    {
        edge = 100;
    }

    else if ( edge >= 2000 )
    {
        edge = 2000;
    }

    do
    {
        u = RandomDouble(-edge,edge);
        v = RandomDouble(0,cp);

        x = cp*exp(-cp_prime*pow(fabs(u),forder));
    }
    while ( v > x );

    return (mean + stddev * u);
}

/*
   Return random integer within a range, lower -> upper INCLUSIVE
*/
int RandomInt(int lower,int upper)
{
   return ((int)(RandomUniform() * (upper - lower + 1)) + lower);
}

/*
   Return random float within a range, lower -> upper
*/
double RandomDouble(double lower,double upper)
{
   return 2*(((upper/2) - (lower/2)) * RandomUniform() + (lower/2));
}


