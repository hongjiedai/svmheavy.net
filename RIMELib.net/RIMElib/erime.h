
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

#ifndef _erime_h
#define _erime_h

/*********************************************************************

RuntIme Mathematical Equation Library
=====================================

Written by: Alistair Shilton
            Melbourne University
            apsh@ecr.mu.oz.au

Date of most recent changes: 13/3/2008

The following is intended for parsing, evaluating, differentiating and
integrating arbitrary (text based) equations.  In normal operation, a
text copy of the equation is given to parseMaths, which parses the
string and returns a pointer to a structure representing that
equation.

Using this structure, one can then evaluate the equation for given
values, differentiate it or perform other operations.  To obtain a
text version of an equation, the function deparseMaths may be called.


Equation string format
======================

Equation strings must be structured in a particular way.  The basic
form of any equation _1_*FnName(_2_) (although there are shortcuts)
is:

FnName[_1_](_2_)

where: - _1_ is a number.  This must be in the usual c standard
         format, or one of the strings:
         inf          = positive infinity
         -inf         = negative infinity
         FIN_INDET    = finite number, exact value unknown
         -FIN_INDET   =  "
         INFIN_INDET  = finite or infinite number, exact value unknown
         -INFIN_INDET =  "
         t_           = (where _ is an integer) error code _

       - _2_ is a list of arguments (possibly empty) of the same form
         FnName[_3_](_4_), separated by commas.

So, using this format, 2*(sin(...)) would be sin[2](...).  Because
this is cumbersome, there are shortcuts.  Firstly, both the [_1_] and
(_2_) arguments are optional (although the latter will be required
unless the function takes no arguments - eg. pi is allowed on its own,
cos is not).

Furthermore the unary operators +,- and the binary operators +,-,*,/,^
are allowed and have the usual meanings and orders of precedence (ie.
posation and negation first, ^ next right to left, the *,/ with equal
importance from left to right, and finally +,- from left to right.  Of
course, brackets have highest precedence).

Finally, a number will be recognised as such, so there is no need to
write one[2.11] - 2.11 will work just fine.

So, the following is a perfectly valid equation:

-sin(pi/2)

Variables are elements of a 2*2 array.  Thus, var(i,j) means "get
variable from position i,j in array given at evaluation time.  So, if
you choose var(1,1)=x then sin(x^2)-acosech(pi+x) must be written:

sin(var(1,1)^2)-acosech(pi+var(1,1))

Alternatively, the shortcuts x = var(1,1), y = var(1,2) and
z = var(1,3) may be used.  A complete list of available functions is
given below.


Local Variables, Placeholders and substitution markers.
=======================================================

var(i,j) actually only insists that i and j are non-negative.  If i
or j is zero (or both) then the variable has a special meaning.  Note
that derivatives are only possibly w.r.t. i>=1 and j>=1.  Also, only
variables i>=1 and j>=1 can be set externally.

For sums, products and integrals, there are local variables defined.
These cannot be set externally, but may be set using the sum, prod
and integ functions, and subsequently used inside said.  They take
the form var(0,i), where i >= 1.

For temporary placeholders, var(i,0) where i>=1 may be used.  If
evaluated, this will give NaN.  However, it is not used internally,
and may be a useful placeholder for later substitutions.

Finally, var(0,0) is used as a placeholder for the substitution
functions, as will be described.  If evaluated it will give NaN.

Note: *DO NOT* reuse local variables recursively - this will break
      badly.  As a rule, this code ASSUMES that recursive re-use does
      not occur, and no guarentees can be made if this rule is
      ignored.  For eg, the following recursively reuses var(0,1)
      and should not be used:
      sum(1,1,10,2,var(0,1)*sum(1,1,10,2,var(0,1)))


Functions: parsing, deparsing, copying and deleting
===================================================

- MathsNode *parseMaths(char *express)

  This function converts the string express (which is assumed to be
  an equation in the correct format, as given above) into a structure
  that may be used by other functions, and return a pointer to this
  structure.  On failure, it will return NULL.

- MathsNode *getZero(void)

  Get "0" quickly.  Will return NULL on failure.

- char *deparseMaths(MathsNode *equation)

  The reverse of parseMaths, this function converts a structural
  representation of an equation into a string (format described above)
  and returns a pointer to this string.  Will return a NULL on
  failure.

- char *deparseMaths_cstruct(MathsNode *equation, char *name)

  Sometimes, it may be desirable to have a static structural
  representation of an equation for use in programming rather than
  re-parsing an equation whenever it is needed in code.  This function
  will return a string containing the necessary code to statically
  initialize such a structure - just printf the result and paste into
  code as required.  The result will have the general form:

  ...
  MathsNode *name = ...;

  where ... is filled appropriately, and name is the string given by
  the calling function.  This variable will point to a static version
  of the equation.  Note that variables of the form MathsNode
  name_content_*_, where * is an integer, and MathsNode *name_refer_*_
  will also be defined.  These are the nodes of the static structure
  and should not be touched.

  Note that the equation so declared must be treated with care.  It
  cannot be modified in any way (eg. differentiation, substitution
  etc) or freed, as the memory is on the stack (or in pre-malloced
  memory if static/global) and hence cannot be dealloced (which these
  functions will attempt).  Exactly what will happen if any of this
  is attempted is architecture dependant, but it's unlikely to be
  pretty.

  This function will return NULL on failure.

- MathsNode *copyEquation(MathsNode *source)

  Create an exact duplicate of the source structure.  Will return a
  NULL on failure.

- int delEqn(MathsNode *equation)

  Free memory associated with equation structure.  Will return 0 on
  success, 1 on failure.  Currently, this operation cannot fail.


Functions: manipulation
=======================

- int simplifyEqn(MathsNode *equation)

  Simplify the equation.  This function is not necessary unless a fast
  function has been called, as simplifyMaths is called after all
  relevant standard functions.  Note that the key assumption made
  during this operation is that 0*f(x) = 0, where x and f are both
  unknown.  This may not be valid.  Will return 0 on success or 1 on
  failure.

- int subVar(MathsNode *equation, MathsNode *var_replace, long i, long j, int naive)

  Replace all instances of var(i,j) in equation with var_replace. If
  naive is nonzero then it will only replace if i and j are known
  integers.  ie. var(5,4) is fair game, var(5,var(8,7)) is not.  On
  success, will return 0.  On failure, will return 1.

- int subVarRow(MathsNode *equation, MathsNode *var_replace, long i, int naive)

  Replace all instances of var(i,~) in equation with var_replace,
  where var(0,0) in var_replace is replaced by whatever ~ is in any
  given instance.  The setting naive in this case acts only on the
  first argumant, i. On success, will return 0.  On failure, will
  return 1.

- int subVarCol(MathsNode *equation, MathsNode *var_replace, long j, int naive)

  Replace all instances of var(~,j) in equation with var_replace,
  where var(0,0) in var_replace is replaced by whatever ~ is in any
  given instance.  The setting naive in this case acts only on the
  second argumant, j. On success, will return 0.  On failure, will
  return 1.

"Fast" versions of the substitution functions are also available.  The
difference between these and standard versions of the functions is
that they do not attempt to simplify the resulting MathsNode in any
way (standard versions do).  This is useful if a lot of calls to these
functions are made sequentially, as in this situation it may be
quicker to simplify things all at once rather than incrementally.  The
fast functions are:

int fast_subVar(MathsNode *equation, MathsNode *var_replace, long i, long j, int naive)
int fast_subVarRow(MathsNode *equation, MathsNode *var_replace, long i, int naive)
int fast_subVarCol(MathsNode *equation, MathsNode *var_replace, long j, int naive)



Functions: differentiation and integration
==========================================

- MathsNode *makeDeriv(MathsNode *equation, long i, long j)

  Compute the partial derivative of equation with respect to var(i,j).
  Will return a pointer to a new structure representing this deriv
  upon success or NULL upon failure.

- MathsNode *MatrixMakeDeriv(MathsNode *equation, long i, long j, long k, long l)

  This form is like the above, except that in this case the derivative
  is wrt var(var(i,j),var(k,l)).  Thus this is a "general" derivative,
  whereas the previous form is a "specific" derivative.  Will return
  NULL upon failure.

- MathsNode *RowMakeDeriv(MathsNode *equation, long i, long j, long k)

  This form is like the above, except that in this case the derivative
  is wrt var(i,var(j,k)).  Will return NULL upon failure.

- MathsNode *ColMakeDeriv(MathsNode *equation, long i, long j, long k)

  This form is like the above, except that in this case the derivative
  is wrt var(var(i,j),k).  Will return NULL upon failure.

Definite integration is also possible using the functions "integ" using
substitution functions.

"Fast" versions of the substitution functions are also available.  The
difference between these and standard versions of the functions is
that they do not attempt to simplify the resulting MathsNode in any
way (standard versions do).  This is useful if a lot of calls to these
functions are made sequentially, as in this situation it may be
quicker to simplify things all at once rather than incrementally.  The
fast functions are:

MathsNode *fast_makeDeriv(MathsNode *equation, long i, long j)
MathsNode *fast_MatrixMakeDeriv(MathsNode *equation, long i, long j, long k, long l)
MathsNode *fast_RowMakeDeriv(MathsNode *equation, long i, long j, long k)
MathsNode *fast_ColMakeDeriv(MathsNode *equation, long i, long j, long k)

  

Functions: evaluation
=====================

The following functions will evaluate an equation.  The value assigned
to var(i,j) is dependent on the function used.  In all cases, if the
result is an error then NaN will be returned.  inf and -inf are also
valid return values here (and valid inputs, too).

- double evaluateEqnMatrix(MathsNode *equation, double **GlobInput)

  Evaluate equation with var(i,j) = GlobInput[i-1][j-1].  No bounds
  checking will be done, so be careful.

- double evaluateEqnVector(MathsNode *equation, double *GlobInput)

  Evaluate equation with var(i,j) = GlobInput[j-1] if i == 1, NaN
  otherwise.  Again, no bounds checking will be done.

- double evaluateEqnList(MathsNode *equation, int numargs, ...)

  This is similar to evaluateEqnVector, except that the vector is
  constructed from the arguments ..., or which there must be numargs.
  Also, var(1,j) = NaN if j > numargs.  However, numargs must be the
  number of arguments in the list ....  Arguments in the list ... must
  be doubles.

- double evaluateEqnGeneral(MathsNode *equation, double (*DoubleArgEvaluationFunction)(void *argContents, long i, long j), void *GlobInput)

  This is the most general form of evaluation.  var(i,j) is retrieved
  by passing i,j and GlobInput (which is assumed to point to some data
  structure that can be understood by DoubleArgEvaluationFunction) to
  the function DoubleArgEvaluationFunction.  This user provided
  function uses these arguments to produce a result x, which is then
  used in var(i,j) = x.

- double evaluateEqnNull(MathsNode *equation)

  Evaluate equation with var(i,j) = NaN.  So, unless the equation is
  just a number, this will trivially give NaN.  It is included for
  completeness.

There are also _e versions of all the above functions.  These are
essentially the same, except that they use call be reference and return
an error code if anything illegal happens (the above just coerce the result
to NaN in this case).

- int evaluateEqnNull_e(double *result, MathsNode *equation)
- int evaluateEqnMatrix_e(double *result, MathsNode *equation, double **GlobInput)
- int evaluateEqnVector_e(double *result, MathsNode *equation, double *GlobInput)
- int evaluateEqnList_e(double *result, MathsNode *equation, int numargs, ...)
- int evaluateEqnGeneral_e(double *result, struct mathsNode *equation, double (*DoubleArgEvaluationFunction)(void *argContents, long i, long j), void *GlobInput)



Functions: hacking
==================

If you add new mathematical functions to this library or make changes
that affect the static data otherwise (eg. correct a derivative) then
the following function may be called to re-calculate the static data
and dump it.  Just pipe the output to a text file, check it, and copy
it over the existing file erime.dat.  Be sure to define INIT_EMATHS
when doing this.

- int initMaths(void)

  This function recompiles all static data (generates erime.dat).


Supported Mathematical Functions
================================

The following standard naming of variables will be used in this
section:

i,j,k,l = integers (no derivative taken w.r.t these).
a,b,c,d = reals (no derivative taken w.r.t these).
x,y,z   = variables.
(finite) means arguments must be finite

Available functions are:

CONSTANTS (0 arguments): one         = 1
                         zero        = 0
                         inf         = inf
                         ninf        = -inf
                         FiniteIndet = FIN_INDET
                         infintIndet = INFIN_INDET
                         isError     = error
VARIABLES (2 arguments): var(i,j)    = variable from array position i,j
          (0 arguments): x           = variable from array position 1,1
                         y           = variable from array position 1,2
                         z           = variable from array position 1,3

NUMBERS (2 arguments): randn(a,b)   = uniformly distr number b/w a and b
                       irand(i,j)   = uniformly distr random int b/w i and j
                       grand(a,b)   = gaussian rand # mean a and stddev b
        (3 arguments): prand(a,b,i) = poly deg i rand # mean a and stddev b

FUNCTIONS:

split(a,b,x,y)  :- return x if a <= b, y otherwise
esplit(a,b,x,y) :- return x if a == b, y otherwise

isint(i,x,y)     :- return x if i = integer, y otherwise
isreal(a,x,y)    :- return x if a = real \ integer, y otherwise
isanum(a,x,y)    :- return x if a = real, y otherwise
isfinint(i,x,y)  :- return x if i = integer (finite), y otherwise
isfinreal(a,x,y) :- return x if a = real \ integer (finite), y otherwise
isfinanum(a,x,y) :- return x if a = real (finite), y otherwise

max(x,y) :- max(x,y) = split(x,y,y,x)
min(x,y) :- min(x,y) = split(x,y,x,y)

inv(x) :- multiplicative inverse of x
neg(x) :- additive inverse of x

add(x,y)  :- x+y
mul(x,y)  :- x*y
div(x,y)  :- x/y
idiv(i,j) :- integer division i/j
fmod(i,j) :- i%j

eval(l,a,x)      :- evaluate x with var(0,l) = a
sum(l,a,b,c,x)   :- sum x over a <= var(0,l) <= b (finite), step size c
prod(l,a,b,c,x)  :- multiply x over a <= var(0,l) <= b (finite), step size c
integ(l,a,b,c,x) :- integrates x w.r.t. var(0,l) from a to b (finite)
                    c > 0 (finite) is the scaled increment size if exact
                    integral not possible (# steps in 1/c).  If c == 0 then
                    DEFAULT_INTEGRAL_STEP will be used instead.

kronDelta(i,j)   :- return 1 if i == j (finite), 0 otherwise
diracDelta(x,y)  :- return posInfty if x ~ y (finite), 0 otherwise
finDiracDel(x,y) :- return 1 if x ~ y (finite), 0 otherwise

abs(x)   :- magnitude of x
sgn(x)   :- sign of x
rint(x)  :- nearest integer (default to higher if borderline)
ceil(x)  :- nearest integer less than
floor(x) :- nearest integer greater than

perm(i,j) :- i perm j = i!/j!
comb(i,j) :- i choose j = i!/(j! (i-j)!)
fact(i)   :- i!

pow(x,y) :- x^y
sqrt(x)  :- square root

sin(x),  cos(x),  tan(x),  cosec(x),  sec(x),  cot(x)
asin(x), acos(x), atan(x), acosec(x), asec(x), acot(x)

sinh(x),  cosh(x),  tanh(x),  cosech(x),  sech(x),  coth(x)
asinh(x), acosh(x), atanh(x), acosech(x), asech(x), acoth(x)

gd(x)  :- the gudermannian: gd(x) = 2atan(exp(x))-pi/2
agd(x) :- the inverse gudermannian: log(tan((x/2)+(pi/4)))

exp(x)   :- exponential
tenup(x) :- 10^x
log(x)   :- natural logarithm log(x)
logX(x)  :- logarithm base 10 log10(x)

gamma(x)     :- gamma function
lngamma(x)   :- log(gamma(x))
gami(y,x)    :- incomplete gamma fn = integ(0 to x) (e^(y-1) * e^-x) dt
gamic(y,x)   :- inverse "  gamma fn = integ(x to inf) (e^(y-1) * e^-x) dt
psi(x)       :- Digamma function = (deriv gamma)(x) / gamma(x)
psi_n(i,x)   :- Polygamma function = ith derivative of psi function.

erf(x)         :- error function
erfc(x)        :- complementary error function
dawson(x)      :- Dawson's fn: D(x) = e^(-x^2) * integ(0 to x)(e^(t^2))dt
normDistr(x)   :- Normal distribution (0 mean, unit variance)
polyDistr(i,x) :- Polynomial distribution (0 mean, unit variance, degree i)


*********************************************************************/

#include <stdarg.h>
#include "gvars.h"
#include "mathtext.h"

#define DEFAULT_INTEGRAL_STEP   1e-4
#define MAX_NUM_LOC_VARS        1024

typedef struct mathsNode
{
    int CombinerFnNumber;
    GenVar multiplier;
    struct mathsNode **SubFunctions;
}
MathsNode;

MathsNode *parseMaths(char *express);
MathsNode *getZero(void);
char *deparseMaths(MathsNode *equation);
char *deparseMaths_cstruct(MathsNode *equation, char *name);
MathsNode *copyEquation(MathsNode *source);
int delEqn(MathsNode *equation);
MathsNode *fast_parseMaths(char *express);

double evaluateEqnNull(MathsNode *equation);
double evaluateEqnMatrix(MathsNode *equation, double **GlobInput);
double evaluateEqnVector(MathsNode *equation, double *GlobInput);
double evaluateEqnList(MathsNode *equation, int numargs, ...);
double evaluateEqnGeneral(struct mathsNode *equation, double (*DoubleArgEvaluationFunction)(void *argContents, long i, long j), void *GlobInput);

int evaluateEqnNull_e(double *result, MathsNode *equation);
int evaluateEqnMatrix_e(double *result, MathsNode *equation, double **GlobInput);
int evaluateEqnVector_e(double *result, MathsNode *equation, double *GlobInput);
int evaluateEqnList_e(double *result, MathsNode *equation, int numargs, ...);
int evaluateEqnGeneral_e(double *result, struct mathsNode *equation, double (*DoubleArgEvaluationFunction)(void *argContents, long i, long j), void *GlobInput);

int simplifyEqn(MathsNode *equation);
int subVar(MathsNode *equation, MathsNode *var_replace, long i, long j, int naive);
int subVarRow(MathsNode *equation, MathsNode *var_replace, long i, int naive);
int subVarCol(MathsNode *equation, MathsNode *var_replace, long j, int naive);
int fast_subVar(MathsNode *equation, MathsNode *var_replace, long i, long j, int naive);
int fast_subVarRow(MathsNode *equation, MathsNode *var_replace, long i, int naive);
int fast_subVarCol(MathsNode *equation, MathsNode *var_replace, long j, int naive);

MathsNode *makeDeriv(MathsNode *start, long i, long j);
MathsNode *MatrixMakeDeriv(MathsNode *start, long i, long j, long k, long l);
MathsNode *RowMakeDeriv(MathsNode *start, long i, long j, long k);
MathsNode *ColMakeDeriv(MathsNode *start, long i, long j, long k);
MathsNode *fast_makeDeriv(MathsNode *equation, long i, long j);
MathsNode *fast_MatrixMakeDeriv(MathsNode *equation, long i, long j, long k, long l);
MathsNode *fast_RowMakeDeriv(MathsNode *equation, long i, long j, long k);
MathsNode *fast_ColMakeDeriv(MathsNode *equation, long i, long j, long k);

int initMaths(void);



#define EERR_BADARGS    MERR_MAX_ERROR-1
#define EERR_MALLOCFAIL MERR_MAX_ERROR-2
#define EERR_INVALIDVAR MERR_MAX_ERROR-3
#define EERR_BADCOPY    MERR_MAX_ERROR-4

#define EERR_MAX_ERROR  MERR_MAX_ERROR-10

#define ERIME_NAN       GSL_NAN

#endif
