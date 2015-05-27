
/*
 *  RIMElib: RuntIme Mathematical Equation Library (underlying functions)
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


#ifndef gmaths_h
#define gmaths_h

#include <math.h>
#include "gvars.h"

GenVar Gen_kronDelta(GenVar ii, GenVar jj);
GenVar Gen_diracDelta(GenVar xx, GenVar yy);
GenVar Gen_finDiracDelta(GenVar xx, GenVar yy);

GenVar Gen_max(GenVar xx, GenVar yy);
GenVar Gen_min(GenVar xx, GenVar yy);

GenVar Gen_abs(GenVar xx);
GenVar Gen_sgn(GenVar xx);
GenVar Gen_rint(GenVar xx);
GenVar Gen_ceil(GenVar xx);
GenVar Gen_floor(GenVar xx);

GenVar Gen_perm(GenVar ii, GenVar jj);
GenVar Gen_comb(GenVar ii, GenVar jj);
GenVar Gen_fact(GenVar ii);

GenVar Gen_pow(GenVar xx, GenVar yy);
GenVar Gen_sqrt(GenVar xx);

GenVar Gen_sin(GenVar xx);
GenVar Gen_cos(GenVar xx);
GenVar Gen_tan(GenVar xx);
GenVar Gen_cosec(GenVar xx);
GenVar Gen_sec(GenVar xx);
GenVar Gen_cot(GenVar xx);
GenVar Gen_asin(GenVar xx);
GenVar Gen_acos(GenVar xx);
GenVar Gen_atan(GenVar xx);
GenVar Gen_acosec(GenVar xx);
GenVar Gen_asec(GenVar xx);
GenVar Gen_acot(GenVar xx);

GenVar Gen_sinh(GenVar xx);
GenVar Gen_cosh(GenVar xx);
GenVar Gen_tanh(GenVar xx);
GenVar Gen_cosech(GenVar xx);
GenVar Gen_sech(GenVar xx);
GenVar Gen_coth(GenVar xx);
GenVar Gen_asinh(GenVar xx);
GenVar Gen_acosh(GenVar xx);
GenVar Gen_atanh(GenVar xx);
GenVar Gen_acosech(GenVar xx);
GenVar Gen_asech(GenVar xx);
GenVar Gen_acoth(GenVar xx);

GenVar Gen_sinc(GenVar xx);
GenVar Gen_gd(GenVar xx);
GenVar Gen_agd(GenVar xx);

GenVar Gen_exp(GenVar xx);
GenVar Gen_tenup(GenVar xx);
GenVar Gen_log(GenVar xx);
GenVar Gen_logX(GenVar xx);

GenVar Gen_gamma(GenVar xx);
GenVar Gen_lngamma(GenVar xx);
GenVar Gen_psi(GenVar xx);
GenVar Gen_psi_n(GenVar ii, GenVar xx);
GenVar Gen_gami(GenVar yy, GenVar xx);
GenVar Gen_gamic(GenVar yy, GenVar xx);

GenVar Gen_erf(GenVar xx);
GenVar Gen_erfc(GenVar xx);
GenVar Gen_normDistr(GenVar xx);
GenVar Gen_polyDistr(GenVar nn, GenVar xx);

GenVar Gen_dawson(GenVar xx);

/*

Gen_kronDelta(i,j)     - standard kronecker delta function
Gen_diracDelta(x,y)    - standard dirac delta function of x-y
Gen_findiracDelta(x,y) - dirac delta of x-y, but 1 rather than inf if x==y

Gen_max(x,y) - return the larger of x and y
Gen_min(x,y) - return the smaller of x and y

Gen_abs(x)   - |x|
Gen_sgn(x)   - sgn(x)
Gen_rint(x)  - round to nearest integer
Gen_ceil(x)  - round to closest integer i > x
Gen_floor(x) - round to closest integer i < x

Gen_perm(i,j) - number of permutations of j objects selectable from i
Gen_comb(i,j) - number of combinations of j objects selectable from i
Gen_fact(i)   - i!

Gen_pow(x,y) - x^y
Gen_sqrt(x)  - square root of x

Gen_sin(x)     - sin(x)
Gen_cos(x)     - cos(x)
Gen_tan(x)     - tan(x)
Gen_cosec(x)   - cosec(x)
Gen_sec(x)     - sec(x)
Gen_cot(x)     - cot(x)
Gen_asin(x)    - arcsin(x)
Gen_acos(x)    - arccos(x)
Gen_atan(x)    - arctin(x)
Gen_acosec(x)  - arccosec(x)
Gen_asec(x)    - arcsec(x)
Gen_acot(x)    - arccot(x)
Gen_sinh(x)    - sinh(x)
Gen_cosh(x)    - cosh(x)
Gen_tanh(x)    - tanh(x)
Gen_cosech(x)  - cosech(x)
Gen_sech(x)    - sech(x)
Gen_coth(x)    - coth(x)
Gen_asinh(x)   - arcsinh(x)
Gen_acosh(x)   - arccosh(x)
Gen_atanh(x)   - arctanh(x)
Gen_acosech(x) - arccosech(x)
Gen_asech(x)   - arcsech(x)
Gen_acoth(x)   - arccoth(x)

Gen_sinc(x) - the sinc function - sinc(x) = sin(x)/x
Gen_gd(x)   - the gudermanian - gd(x) = 2*arctan(exp(x)) - pi/2
Gen_agd(x)  - the inverse gudermanian

Gen_exp(x)   - e^x
Gen_tenup(x) - 10^x
Gen_log(x)   - natural logarithm
Gen_logX(x)  - log base 10

Gen_gamma(a)     - gamma function = integ(0 to inf) (e^(a-1) * e^-x) dt
Gen_lngamma(a)   - log gamma function = log(|gamma(x)|)
Gen_psi(x)       - Digamma function = (deriv gamma)(x) / gamma(x)
Gen_psi_n(x)     - Polygamma function = ith derivative of psi(x)
Gen_gami(a,x)    - incomplete gamma fn = integ(0 to x) (e^(a-1) * e^-x) dt
Gen_gamic(a,x)   - inverse "  gamma fn = integ(x to inf) (e^(a-1) * e^-x) dt
Gen_beta(a,b)    - beta function = gamma(a)*gamma(b)/gamma(a+b)
Gen_lbeta(a,b)   - log beta function = log(beta(a,b))
Gen_betaI(x,a,b) - incompl bet fn = integ(0 to x) (t^(a-1) * (1-t)^(b-1)) dt

Gen_erf(x)         - error function
Gen_erfc(x)        - complementary error function
Gen_dawson(x)      - Dawson's fn: D(x) = e^(-x^2) * integ(0 to x)(e^(-t^2))dt
Gen_normDistr(x)   - Normal distribution (0 mean, unit variance)
Gen_polyDistr(n,x) - Polynomial distribution (0 mean, unit variance)

*/

#endif
