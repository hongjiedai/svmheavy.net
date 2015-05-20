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

#include <iostream>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "gsl/gsl_math.h"

#include "common/svdefs.h"
#include "common/vector.h"
#include "common/sparsevector.h"
#include "common/matrix.h"
#include "common/friends.h"
#include "common/search.h"
#include "kernel.h"
#include "common/outfilt.h"


#define VAR_REP_FORM_X_GENERAL  "var(5,7)*sum(1,1,var(0,0),1,var(9+var(0,0),var(0,1))*(var(6,var(0,1))-var(8,var(0,1))))"
#define VAR_REP_FORM_Y_GENERAL  "var(5,7)*sum(1,1,var(0,0),1,var(9+var(0,0),var(0,1))*(var(7,var(0,1))-var(8,var(0,1))))"

#define VAR_REP_FORM_X_DIAG     "var(5,7)*var(9+var(0,0),var(0,0))*(var(6,var(0,0))-var(8,var(0,0)))"
#define VAR_REP_FORM_Y_DIAG     "var(5,7)*var(9+var(0,0),var(0,0))*(var(7,var(0,0))-var(8,var(0,0)))"

#define DOT_PRODUCT_NONSPARSE           "sum(11,1,var(5,8),1,var(1,var(0,11))*var(2,var(0,11)))"
#define DIFF_MAG_SQUARED_NONSPARSE      "sum(11,1,var(5,8),1,(var(1,var(0,11))-var(2,var(0,11)))*(var(1,var(0,11))-var(2,var(0,11))))"
#define DIFF_MAGNITUDE_NONSPARSE        "sqrt(sum(11,1,var(5,8),1,(var(1,var(0,11))-var(2,var(0,11)))*(var(1,var(0,11))-var(2,var(0,11)))))"

#define DOT_PRODUCT_SPARSE              "sum(11,2,2*var(6,1),2,sum(12,2,2*var(7,1),2,esplit(var(6,var(0,11)),var(7,var(0,12)),var(6,1+var(0,11))*var(7,1+var(0,12)),0)))"
#define DIFF_MAG_SQUARED_SPARSE         "sum(11,2,2*var(6,1),2,var(6,1+var(0,11))*var(6,1+var(0,11)))+sum(12,2,2*var(7,1),2,var(7,1+var(0,12))*var(7,1+var(0,12)))-(2*sum(11,2,2*var(6,1),2,sum(12,2,2*var(7,1),2,esplit(var(6,var(0,11)),var(7,var(0,12)),var(6,1+var(0,11))*var(7,1+var(0,12)),0))))"
#define DIFF_MAGNITUDE_SPARSE           "sqrt(sum(11,2,2*var(6,1),2,var(6,1+var(0,11))*var(6,1+var(0,11)))+sum(12,2,2*var(7,1),2,var(7,1+var(0,12))*var(7,1+var(0,12)))-(2*sum(11,2,2*var(6,1),2,sum(12,2,2*var(7,1),2,esplit(var(6,var(0,11)),var(7,var(0,12)),var(6,1+var(0,11))*var(7,1+var(0,12)),0)))))"


#define N_KERNS         51
#define MAX_DESCR_LEN   400
#define MAX_DEFN_LEN    400


//
// Descriptions of kernel functions
//

char kern_descr[N_KERNS][MAX_DESCR_LEN] =
{
" 1. Dot product (inbuilt): (x'y)\n",
" 2. Difference measure (inbuilt): (||x-y||^2)\n",
" 3. Euclidean distance (inbuilt): ||x-y||\n",
" 4. Polynomial 1a (inbuilt): (x'y)^uc1\n",
" 5. Polynomial 1b (inbuilt): (x'y)^uv1\n",
" 6. Polynomial 2a (inbuilt): ((x'y)^uc1 + uc2)^uc3\n",
" 7. Polynomial 2a (inbuilt): ((x'y)^uv1 + uc1)^uc2\n",
" 8. Polynomial 2a (inbuilt): ((x'y)^uc1 + uv1)^uc2\n",
" 9. Polynomial 2a (inbuilt): ((x'y)^uc1 + uc2)^uv1\n",
"10. Polynomial 2a (inbuilt): ((x'y)^uc1 + uv1)^uv2\n",
"11. Polynomial 2a (inbuilt): ((x'y)^uv1 + uc1)^uv2\n",
"12. Polynomial 2a (inbuilt): ((x'y)^uv1 + uv2)^uc1\n",
"13. Polynomial 2a (inbuilt): ((x'y)^uv1 + uv2)^uv3\n",
"14. neural network 1 (inbuilt): tanh((x'y))\n",
"15. neural network 2 (inbuilt): tanh((x'y) + uc1)\n",
"16. neural network 2 (inbuilt): tanh((x'y) + uv1)\n",
"17. exponential expansion 1 (inbuilt): exp((x'y))\n",
"18. exponential expansion 2 (inbuilt): exp((x'y)) - uc1\n",
"19. exponential expansion 2 (inbuilt): exp((x'y)) - uv1\n",
"20. gaussian RBF 1 (inbuilt): exp(-(||x-y||^2))\n",
"21. gaussian RBF 2 (inbuilt): exp(-(||x-y||^2)) - uc1\n",
"22. gaussian RBF 2 (inbuilt): exp(-(||x-y||^2)) - uv1\n",
"23. exponential RBF 1 (inbuilt): exp(-||x-y||)\n",
"24. exponential RBF 2 (inbuilt): exp(-||x-y||) - uc1\n",
"25. exponential RBF 2 (inbuilt): exp(-||x-y||) - uv1\n",
"26. sigmoidal (nonsparse): product(i=1 to dim)(1 / ( 1 + exp(-z_i) ))\n",
"27. error function (nonsparse): prod(i=1 to dim)( (1/2) * (1 + erf(z_i) ))\n",
"28. vovk's real polynomial: (1-((x'y)^uc1))/(1-x'y) if x'y != 1, uv1 otherwise\n",
"29. vovk's real polynomial: (1-((x'y)^uv1))/(1-x'y) if x'y != 1, uv1 otherwise\n",
"30. dirichlet (nonsparse): prod(i=1-dim)(sin((uc1+0.5)*z_i)/(2*sin(0.5*z_i)))\n",
"31. dirichlet (nonsparse): prod(i=1-dim)(sin((uv1+0.5)*z_i)/(2*sin(0.5*z_i)))\n",
"32. weak_fourier: pi*cosh(pi-||x-y||)\n",
"33. strong_four (nonsparse): prod(i=1,dim)((1-(uc1^2))/2*((1-(2*uc1*cos(z_i))+\n    (uc1^2))))\n",
"34. strong_four (nonsparse): prod(i=1,dim)((1-(uv1^2))/2*((1-(2*uv1*cos(z_i))+\n    (uv1^2))))\n",
"35. Trig poly degree n (nonsparse): product(i=1 to dim)(sin(uc1+(1/2)).\n    (xm_i-ym_i)/sin((xm_i-ym_i)/2))\n",
"36. Trig poly degree n (nonsparse): product(i=1 to dim)(sin(uv1+(1/2)).\n    (xm_i-ym_i)/sin((xm_i-ym_i)/2))\n",
"37. Lin spline with infinite # points (nonsparse): prod(i=1-dim)(1+xm_i*ym_i +\n    xm_i*ym_i*min(xm_i,ym_i) - 0.5*(xm_i+ym_i)*((min(xm_i,ym_i))^2) +\n    ((min(xm_i,ym_i))^3)/3)\n",
"38. dense_spline (non''): prod(i=1-dim)((1/2).min(xm_i,ym_i).(max(xm_i,ym_i)^2)\n    - (1/6).(max(xm_i,ym_i)^3) + xm_i.ym_i)\n",
"39. const_spline (nonsparse): product(i=1 to dim)(k_i)\n             k_i = 0                         , if xm_i-ym_i < -uc1\n                   ((xm_i-ym_i)+uc1)/(2*uc1) , if -uc1 <= xm_i-ym_i <= uc1\n                   1                         , if xm_i-ym_i > uc1\n",
"40. const_spline (nonsparse): product(i=1 to dim)(k_i)\n             k_i = 0                         , if xm_i-ym_i < -uv1\n                   ((xm_i-ym_i)+uv1)/(2*uv1) , if -uv1 <= xm_i-ym_i <= uv1\n                   1                         , if xm_i-ym_i > uv1\n",
"41. Thin plate Spline 1: (||x-y||^2)^(uc1+(1/2))\n",
"42. Thin plate Spline 1: (||x-y||^2)^(uv1+(1/2))\n",
"43. Thin plate Spline 2: ((||x-y||^2)^uc1).ln((||x-y||^2)^(1/2))\n",
"44. Thin plate Spline 2: ((||x-y||^2)^uv1).ln((||x-y||^2)^(1/2))\n",
"45. B-Splines (nonsparse): product(i=1 to dim)(xm_i - ym_i)\n",
"46. Diagonal: uc1 if xi == yi, 0 otherwise.\n",
"47. Diagonal: uc1 if xi == yi and xd == +1, 0 otherwise.\n",
"48. Diagonal: uc1 if xi == yi and xd == -1, 0 otherwise.\n",
"49. Constant: uc1 - for use in automatically biased kernels.\n",
"50. dirichletb (nonsparse): sum(i=1 to dim)(sin((uc1+0.5)*z_i)/sin(0.5*z_i))\n",
"51. Tube shrinking: xd.yd.uc1.\n"
};


//
// Variable use of inbuilt kernels - uc, uv
//

int inbuilt_var_use[N_KERNS][2] =
{
{0,0},
{0,0},
{0,0},
{1,0},
{0,1},
{3,0},
{2,1},
{2,1},
{2,1},
{1,2},
{1,2},
{1,2},
{0,3},
{0,0},
{1,0},
{0,1},
{0,0},
{1,0},
{0,1},
{0,0},
{1,0},
{0,1},
{0,0},
{1,0},
{0,1},
{0,0},
{0,0},
{1,0},
{0,1},
{1,0},
{0,1},
{0,0},
{1,0},
{0,1},
{1,0},
{0,1},
{0,0},
{0,0},
{0,1},
{0,1},
{0,1},
{0,0},
{1,0},
{1,0},
{1,0},
{1,0},
{1,0},
{1,0}
};


//
// Mathematically passable versions of kernel functions
//

char inbuilt_kernels[N_KERNS][MAX_DEFN_LEN] =
{
"var(5,12)",
"var(5,13)", 
"var(5,14)",
"var(5,12)^var(4,1)",
"var(5,12)^var(3,1)",
"(var(5,12)^var(4,1)+var(4,2))^var(4,3)",
"(var(5,12)^var(3,1)+var(4,1))^var(4,2)",
"(var(5,12)^var(4,1)+var(3,1))^var(4,2)",
"(var(5,12)^var(4,1)+var(4,2))^var(3,1)",
"(var(5,12)^var(4,1)+var(3,1))^var(3,2)",
"(var(5,12)^var(3,1)+var(4,1))^var(3,2)",
"(var(5,12)^var(3,1)+var(3,2))^var(4,1)",
"(var(5,12)^var(3,1)+var(3,2))^var(3,3)",
"tanh(var(5,12))",
"tanh(var(5,12)+var(4,1))",
"tanh(var(5,12)+var(3,1))",
"exp(var(5,12))",
"exp(var(5,12))-var(4,1)",
"exp(var(5,12))-var(3,1)",
"exp(-var(5,13))",
"exp(-var(5,13))-var(4,1)",
"exp(-var(5,13))-var(3,1)",
"exp(-var(5,14))",
"exp(-var(5,14))-var(4,1)",
"exp(-var(5,14))-var(3,1)",
"prod(11,1,var(5,8),1,inv(1+exp(-var(1,var(0,11))+var(2,var(0,11))))",
"prod(11,1,var(5,8),1,(1+erf((1/sqrt(2))*(var(1,var(0,11))-var(2,var(0,11)))))/2)",
"esplit(1-var(5,12),0,var(4,1),(1-var(5,12)^var(4,1))/(1-var(5,12)))",
"esplit(1-var(5,12),0,var(3,1),(1-var(5,12)^var(3,1))/(1-var(5,12)))",
"prod(11,1,var(5,8),1,esplit(sin(var(1,var(0,11))-var(2,var(0,11))),0,0.5+var(4,1),sin((var(1,var(0,11))-var(2,var(0,11)))*(0.5+var(4,1)))/(2*sin(0.5*(var(1,var(0,11))-var(2,var(0,11)))))))",
"prod(11,1,var(5,8),1,esplit(sin(var(1,var(0,11))-var(2,var(0,11))),0,0.5+var(3,1),sin((var(1,var(0,11))-var(2,var(0,11)))*(0.5+var(3,1)))/(2*sin(0.5*(var(1,var(0,11))-var(2,var(0,11)))))))",
"pi*cosh(pi-var(5,14))",
"prod(11,1,var(5,8),1,(1-var(4,1)^2)*0.5/(1+(-2*var(4,1)*cos(var(1,var(0,11))-var(2,var(0,11))))+(var(4,1)^2)))",
"prod(11,1,var(5,8),1,(1-var(3,1)^2)*0.5/(1+(-2*var(3,1)*cos(var(1,var(0,11))-var(2,var(0,11))))+(var(3,1)^2)))",
"prod(11,1,var(5,8),1,esplit(sin(var(1,var(0,11))-var(2,var(0,11))),0,(sin(var(4,1)+0.5)*(var(1,var(0,11))-var(2,var(0,11)))/(sin(0.5*(var(1,var(0,11))-var(2,var(0,11)))))),(2*sin(var(4,1)+0.5))))",
"prod(11,1,var(5,8),1,esplit(sin(var(1,var(0,11))-var(2,var(0,11))),0,(sin(var(3,1)+0.5)*(var(1,var(0,11))-var(2,var(0,11)))/(sin(0.5*(var(1,var(0,11))-var(2,var(0,11)))))),(2*sin(var(3,1)+0.5))))",
"prod(11,1,var(5,8),1,1+(var(1,var(0,11))*var(2,var(0,11)))+((var(1,var(0,11))*var(2,var(0,11)))*min(var(1,var(0,11)),var(2,var(0,11))))-(0.5*(var(1,var(0,11))+var(2,var(0,11)))*(min(var(1,var(0,11)),var(2,var(0,11)))^2))+((min(var(1,var(0,11)),var(2,var(0,11)))^3)/3))",
"prod(11,1,var(5,8),1,(0.5*min(var(1,var(0,11)),var(2,var(0,11)))*(max(var(1,var(0,11)),var(2,var(0,11)))^2)) - ((1/6)*(max(var(1,var(0,11)),var(2,var(0,11)))^3)) + (var(1,var(0,11))*var(2,var(0,11))))",
"prod(11,1,var(5,8),1,split((var(1,var(0,11))-var(2,var(0,11))),-var(4,1),0,split((var(1,var(0,11))-var(2,var(0,11))),var(4,1),(((var(1,var(0,11))-var(2,var(0,11)))+var(4,1))/(2*var(4,1))),1)))",
"prod(11,1,var(5,8),1,split((var(1,var(0,11))-var(2,var(0,11))),-var(3,1),0,split((var(1,var(0,11))-var(2,var(0,11))),var(3,1),(((var(1,var(0,11))-var(2,var(0,11)))+var(3,1))/(2*var(3,1))),1)))",
"var(5,13)^(var(4,1)+0.5)",
"var(5,13)^(var(3,1)+0.5)",
"(var(5,13)^var(4,1))*log(sqrt(var(5,13)))",
"(var(5,13)^var(3,1))*log(sqrt(var(5,13)))",
"prod(11,1,var(5,8),1,var(1,var(0,11))-mvar(2,var(0,11)))",
"esplit(var(5,3),var(5,4),var(4,1),0)",
"esplit(var(5,3),var(5,4),esplit(var(5,1),1,var(4,1),0),0)",
"esplit(var(5,3),var(5,4),esplit(var(5,1),-1,var(4,1),0),0)",
"var(4,1)",
"sum(11,1,var(5,8),1,esplit(sin(0.5*(var(1,var(0,11))-var(2,var(0,11)))),0,(sin((var(1,var(0,11))-var(2,var(0,11)))*(0.5+var(4,1)))*inv(sin(0.5*(var(1,var(0,11)-var(2,var(0,11)))))),(2*(0.5+var(4,1)))))",
"var(4,1)*var(5,1)*var(5,2)"
};



//
// Table of pointers to optimised kernels (or NULL if none available).
//

typedef struct
{
    double (*kernel)(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
    double (*kernel_scale_prod)(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
    double (*kernel_sparse)(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
    double (*kernel_sparse_opt)(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
}
OptimKernel;

double fast_kernel_01_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_02_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_03_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_04_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_05_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_06_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_07_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_08_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_09_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_10_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_11_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_12_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_13_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_14_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_15_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_16_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_17_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_18_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_19_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_20_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_21_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_22_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_23_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_24_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_25_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);

double fast_kernel_01_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_02_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_03_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_04_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_05_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_06_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_07_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_08_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_09_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_10_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_11_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_12_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_13_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_14_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_15_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_16_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_17_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_18_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_19_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_20_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_21_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_22_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_23_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_24_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_25_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);

double fast_kernel_01_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_02_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_03_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_04_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_05_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_06_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_07_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_08_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_09_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_10_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_11_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_12_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_13_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_14_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_15_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_16_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_17_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_18_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_19_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_20_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_21_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_22_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_23_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_24_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
double fast_kernel_25_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);

double fast_opt_kernel_01_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_02_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_03_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_04_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_05_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_06_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_07_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_08_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_09_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_10_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_11_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_12_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_13_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_14_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_15_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_16_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_17_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_18_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_19_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_20_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_21_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_22_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_23_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_24_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
double fast_opt_kernel_25_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);

OptimKernel inbuilt_opt_kernel[N_KERNS] =
{
{ fast_kernel_01_ , fast_kernel_01_scale_prod , fast_kernel_01_sparse , fast_opt_kernel_01_sparse },
{ fast_kernel_02_ , fast_kernel_02_scale_prod , fast_kernel_02_sparse , fast_opt_kernel_02_sparse },
{ fast_kernel_03_ , fast_kernel_03_scale_prod , fast_kernel_03_sparse , fast_opt_kernel_03_sparse },
{ fast_kernel_04_ , fast_kernel_04_scale_prod , fast_kernel_04_sparse , fast_opt_kernel_04_sparse },
{ fast_kernel_05_ , fast_kernel_05_scale_prod , fast_kernel_05_sparse , fast_opt_kernel_05_sparse },
{ fast_kernel_06_ , fast_kernel_06_scale_prod , fast_kernel_06_sparse , fast_opt_kernel_06_sparse },
{ fast_kernel_07_ , fast_kernel_07_scale_prod , fast_kernel_07_sparse , fast_opt_kernel_07_sparse },
{ fast_kernel_08_ , fast_kernel_08_scale_prod , fast_kernel_08_sparse , fast_opt_kernel_08_sparse },
{ fast_kernel_09_ , fast_kernel_09_scale_prod , fast_kernel_09_sparse , fast_opt_kernel_09_sparse },
{ fast_kernel_10_ , fast_kernel_10_scale_prod , fast_kernel_10_sparse , fast_opt_kernel_10_sparse },
{ fast_kernel_11_ , fast_kernel_11_scale_prod , fast_kernel_11_sparse , fast_opt_kernel_11_sparse },
{ fast_kernel_12_ , fast_kernel_12_scale_prod , fast_kernel_12_sparse , fast_opt_kernel_12_sparse },
{ fast_kernel_13_ , fast_kernel_13_scale_prod , fast_kernel_13_sparse , fast_opt_kernel_13_sparse },
{ fast_kernel_14_ , fast_kernel_14_scale_prod , fast_kernel_14_sparse , fast_opt_kernel_14_sparse },
{ fast_kernel_15_ , fast_kernel_15_scale_prod , fast_kernel_15_sparse , fast_opt_kernel_15_sparse },
{ fast_kernel_16_ , fast_kernel_16_scale_prod , fast_kernel_16_sparse , fast_opt_kernel_16_sparse },
{ fast_kernel_17_ , fast_kernel_17_scale_prod , fast_kernel_17_sparse , fast_opt_kernel_17_sparse },
{ fast_kernel_18_ , fast_kernel_18_scale_prod , fast_kernel_18_sparse , fast_opt_kernel_18_sparse },
{ fast_kernel_19_ , fast_kernel_19_scale_prod , fast_kernel_19_sparse , fast_opt_kernel_19_sparse },
{ fast_kernel_20_ , fast_kernel_20_scale_prod , fast_kernel_20_sparse , fast_opt_kernel_20_sparse },
{ fast_kernel_21_ , fast_kernel_21_scale_prod , fast_kernel_21_sparse , fast_opt_kernel_21_sparse },
{ fast_kernel_22_ , fast_kernel_22_scale_prod , fast_kernel_22_sparse , fast_opt_kernel_22_sparse },
{ fast_kernel_23_ , fast_kernel_23_scale_prod , fast_kernel_23_sparse , fast_opt_kernel_23_sparse },
{ fast_kernel_24_ , fast_kernel_24_scale_prod , fast_kernel_24_sparse , fast_opt_kernel_24_sparse },
{ fast_kernel_25_ , fast_kernel_25_scale_prod , fast_kernel_25_sparse , fast_opt_kernel_25_sparse },
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL},
{NULL,NULL,NULL,NULL}
};



//
// dot_prod: xm'ym
// dif_meas: ||xm-ym||^2 = (xm-ym)'(xm-ym)
// eucldist: ||xm-ym||
//
// Non-sparse: xm = covw.covar.(x-mean)
//             ym = covw.covar.(y-mean)
//
// Sparse: xm = x
//         ym = x
//

double dot_prod(const fVECTOR &x, const fVECTOR &y, double covw, const fVECTOR &mean, const fMATRIX &covar);
double dif_meas(const fVECTOR &x, const fVECTOR &y, double covw, const fVECTOR &mean, const fMATRIX &covar);
double eucldist(const fVECTOR &x, const fVECTOR &y, double covw, const fVECTOR &mean, const fMATRIX &covar);

double dot_prod_scale_prod(const fVECTOR &x, const fVECTOR &y, double covw, const fVECTOR &mean, const fMATRIX &covar);
double dif_meas_scale_prod(const fVECTOR &x, const fVECTOR &y, double covw, const fVECTOR &mean, const fMATRIX &covar);
double eucldist_scale_prod(const fVECTOR &x, const fVECTOR &y, double covw, const fVECTOR &mean, const fMATRIX &covar);




/******************************************************************/
/*                                                                */
/* Simple kernel class                                            */
/*                                                                */
/******************************************************************/

//
// Constructors and Destructors
// ============================
//
// default: uncorrected fixed linear kernel, covw = 1, scale_prod = 0
// copy: copy the source
// complete: all variables must be provided.
//

basic_kernel::basic_kernel()
{
    try
    {
        where_to   = NULL;
        echo_level = 0;

        kernel_def = LINEAR_KERNEL;

        is_cross_kernel_valid = 0;
        is_kernel_deriv_valid = 0;

        scale_prod    = 0;
        sparse_kernel = 0;

        covw = 1;
        covar.make_scalar(DO_FORCE);

        fix_covw  = 1;
        fix_mean  = 1;
        fix_covar = 1;

        is_opt            = 0;
        optim_kernel      = NULL;
        optim_kernel_fast = NULL;
    }

    catch ( int x )
    {
        L_THROW(x);
    }

    simplify_kernel();

    return;
}

basic_kernel::basic_kernel(const basic_kernel &what)
{
    try
    {
        where_to   = NULL;
        echo_level = 0;

        kernel_def       = what.kernel_def;
        cross_kernel_def = what.cross_kernel_def;
        kernel_deriv_def = what.kernel_deriv_def;

        is_cross_kernel_valid = what.is_cross_kernel_valid;
        is_kernel_deriv_valid = what.is_kernel_deriv_valid;

        scale_prod    = what.scale_prod;
        sparse_kernel = what.sparse_kernel;

        uc = what.uc;
        uv = what.uv;

        covw  = what.covw;
        mean  = what.mean;
        covar = what.covar;

        fix_covw  = what.fix_covw;
        fix_mean  = what.fix_mean;
        fix_covar = what.fix_covar;

        is_opt            = what.is_opt;
        optim_kernel      = what.optim_kernel;
        optim_kernel_fast = what.optim_kernel_fast;
    }

    catch ( int x )
    {
        L_THROW(x);
    }

    return;
}

basic_kernel::basic_kernel(char *kern_def, const fVECTOR &_uc, const fVECTOR &_uv, double _covw, const fVECTOR &_mean, const fMATRIX &_covar, int _fix_covw, int _fix_mean, int _fix_covar, int _scale_prod, int _sparse_kernel)
{
    long i,j;

    try
    {
        where_to   = NULL;
        echo_level = 0;

        kernel_def = kern_def;

        is_cross_kernel_valid = 0;
        is_kernel_deriv_valid = 0;

        scale_prod    = _scale_prod;
        sparse_kernel = _sparse_kernel;

        uc = _uc;
        uv = _uv;

        covw  = _covw;
        mean  = _mean;
        covar = _covar;

        fix_covw  = _fix_covw;
        fix_mean  = _fix_mean;
        fix_covar = _fix_covar;

        j = uc.get_effective_size();

        for ( i = 1 ; i <= j ; i++ )
        {
            kernel_def.naive_replace(uc[i-1],4,i);
        }

        is_opt            = 0;
        optim_kernel      = NULL;
        optim_kernel_fast = NULL;
    }

    catch ( int x )
    {
        L_THROW(x);
    }

    simplify_kernel();

    return;
}

basic_kernel::basic_kernel(int kern_number, const fVECTOR &_uc, const fVECTOR &_uv, double _covw, const fVECTOR &_mean, const fMATRIX &_covar, int _fix_covw, int _fix_mean, int _fix_covar, int _scale_prod, int _sparse_kernel)
{
    long i,j;

    THROW_ASSERT(kern_number >= 1);
    THROW_ASSERT(kern_number <= N_KERNS);

    try
    {
        where_to   = NULL;
        echo_level = 0;

        kernel_def = inbuilt_kernels[kern_number-1];

        is_cross_kernel_valid = 0;
        is_kernel_deriv_valid = 0;

        scale_prod    = _scale_prod;
        sparse_kernel = _sparse_kernel;

        uc = _uc;
        uv = _uv;

        covw  = _covw;
        mean  = _mean;
        covar = _covar;

        fix_covw  = _fix_covw;
        fix_mean  = _fix_mean;
        fix_covar = _fix_covar;

        j = uc.get_effective_size();

        for ( i = 1 ; i <= j ; i++ )
        {
            kernel_def.naive_replace(uc[i-1],4,i);
        }

        is_opt            = 0;
        optim_kernel      = NULL;
        optim_kernel_fast = NULL;

        if ( sparse_kernel )
        {
            if ( (inbuilt_opt_kernel[kern_number-1]).kernel_sparse != NULL )
            {
                is_opt            = kern_number;
                optim_kernel      = (inbuilt_opt_kernel[kern_number-1]).kernel_sparse;
                optim_kernel_fast = (inbuilt_opt_kernel[kern_number-1]).kernel_sparse_opt;
            }
        }

        else
        {
            if ( scale_prod )
            {
                if ( (inbuilt_opt_kernel[kern_number-1]).kernel_scale_prod != NULL )
                {
                    is_opt            = kern_number;
                    optim_kernel      = (inbuilt_opt_kernel[kern_number-1]).kernel_scale_prod;
                    optim_kernel_fast = NULL;
                }
            }

            else
            {
                if ( (inbuilt_opt_kernel[kern_number-1]).kernel != NULL )
                {
                    is_opt            = kern_number;
                    optim_kernel      = (inbuilt_opt_kernel[kern_number-1]).kernel;
                    optim_kernel_fast = NULL;
                }
            }
        }
    }

    catch ( int x )
    {
        L_THROW(x);
    }

    simplify_kernel();

    return;
}

basic_kernel &basic_kernel::operator=(const basic_kernel &source)
{
    try
    {
        if ( this != &source )
        {
            uc.make_zero(DO_FORCE);
            uv.make_zero(DO_FORCE);

            covar = 0;
            mean.make_zero(DO_FORCE);
            covar.make_zero(DO_FORCE);

            kernel_def       = source.kernel_def;
            cross_kernel_def = source.cross_kernel_def;
            kernel_deriv_def = source.kernel_deriv_def;

            is_cross_kernel_valid = source.is_cross_kernel_valid;
            is_kernel_deriv_valid = source.is_kernel_deriv_valid;

            scale_prod    = source.scale_prod;
            sparse_kernel = source.sparse_kernel;

            uc = source.uc;
            uv = source.uv;

            covw  = source.covw;
            mean  = source.mean;
            covar = source.covar;

            fix_covw  = source.fix_covw;
            fix_mean  = source.fix_mean;
            fix_covar = source.fix_covar;

            is_opt            = source.is_opt;
            optim_kernel      = source.optim_kernel;
            optim_kernel_fast = source.optim_kernel_fast;
        }
    }

    catch ( int x )
    {
        L_THROW(x);
    }

    return (*this);
}




//
// Kernel evaluation functions
// ===========================
//
// Calculate the relevant result (after calculating the relevant
// derivatives if this is the first call).
//

double basic_kernel::kernel(const fVECTOR &x, const fVECTOR &y, L_DOUBLE *x_fast, L_DOUBLE *y_fast, double xd, double yd, long xi, long yi, double xz, double yz)
{
    void *argumentus[11];
    double confuseus[16];
    double result;

    if ( is_opt )
    {
        if ( ( x_fast == NULL ) || ( y_fast == NULL ) || ( optim_kernel_fast == NULL ) )
        {
            result = optim_kernel(x,y,xd,yd,xi,yi,xz,yz,uc,uv,covw,mean,covar);
        }

        else
        {
            result = optim_kernel_fast(x_fast,y_fast,&(uc[0]),&(uv[0]),&covw);
        }
    }

    else
    {
        argumentus[1-1]  = NULL;
        argumentus[2-1]  = NULL;
        argumentus[3-1]  = (void *) &uv;
        argumentus[4-1]  = (void *) &uc;
        argumentus[5-1]  = (void *) confuseus;
        argumentus[6-1]  = (void *) &x;
        argumentus[7-1]  = (void *) &y;
        argumentus[8-1]  = (void *) &mean;
        argumentus[9-1]  = NULL;
        argumentus[10-1] = (void *) &covar;
        argumentus[11-1] = (void *) &scale_prod;

        confuseus[1-1]  = (double) xd;
        confuseus[2-1]  = (double) yd;
        confuseus[3-1]  = (double) xi;
        confuseus[4-1]  = (double) yi;
        confuseus[5-1]  = (double) xz;
        confuseus[6-1]  = (double) yz;
        confuseus[7-1]  = (double) covw;

        confuseus[10-1] = GSL_NAN;
        confuseus[11-1] = GSL_NAN;
        confuseus[12-1] = GSL_NAN;
        confuseus[13-1] = GSL_NAN;
        confuseus[14-1] = GSL_NAN;

        if ( sparse_kernel )
        {
            if ( x.get_offset_element(2-1) <= y.get_offset_element(2-1) )
            {
                confuseus[8-1]  = x.get_offset_element(2-1);
                confuseus[9-1]  = y.get_offset_element(2-1);
                confuseus[15-1] = x.get_offset_element(2-1);
                confuseus[16-1] = y.get_offset_element(2-1);
            }

            else
            {
                confuseus[8-1]  = y.get_offset_element(2-1);
                confuseus[9-1]  = x.get_offset_element(2-1);
                confuseus[15-1] = x.get_offset_element(2-1);
                confuseus[16-1] = y.get_offset_element(2-1);
            }
        }

        else
        {
            if ( x.get_effective_size() <= y.get_effective_size() )
            {
                confuseus[8-1]  = (double) x.get_effective_size();
                confuseus[9-1]  = (double) y.get_effective_size();
                confuseus[15-1] = (double) x.get_effective_size();
                confuseus[16-1] = (double) y.get_effective_size();
            }

            else
            {
                confuseus[8-1]  = (double) y.get_effective_size();
                confuseus[9-1]  = (double) x.get_effective_size();
                confuseus[15-1] = (double) x.get_effective_size();
                confuseus[16-1] = (double) y.get_effective_size();
            }
        }

        try
        {
            result = kernel_def(&DoubleArgEvaluationFunctionKernel,(void *) argumentus);
        }

        catch ( int x )
        {
            L_THROW(x);
        }
    }

    if ( !gsl_finite(result) )
    {
        L_THROW(0);
    }

    return result;
}

double basic_kernel::cross_kernel(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz)
{
    void *argumentus[11];
    double confuseus[16];
    double result;

    argumentus[1-1]  = NULL;
    argumentus[2-1]  = NULL;
    argumentus[3-1]  = (void *) &uv;
    argumentus[4-1]  = (void *) &uc;
    argumentus[5-1]  = (void *) confuseus;
    argumentus[6-1]  = (void *) &x;
    argumentus[7-1]  = (void *) &y;
    argumentus[8-1]  = (void *) &mean;
    argumentus[9-1]  = NULL;
    argumentus[10-1] = (void *) &covar;
    argumentus[11-1] = (void *) &scale_prod;

    confuseus[1-1]  = (double) xd;
    confuseus[2-1]  = (double) yd;
    confuseus[3-1]  = (double) xi;
    confuseus[4-1]  = (double) yi;
    confuseus[5-1]  = (double) xz;
    confuseus[6-1]  = (double) yz;
    confuseus[7-1]  = (double) covw;

    confuseus[10-1] = GSL_NAN;
    confuseus[11-1] = GSL_NAN;
    confuseus[12-1] = GSL_NAN;
    confuseus[13-1] = GSL_NAN;
    confuseus[14-1] = GSL_NAN;

    if ( sparse_kernel )
    {
        if ( x.get_offset_element(2-1) <= y.get_offset_element(2-1) )
        {
            confuseus[8-1]  = x.get_offset_element(2-1);
            confuseus[9-1]  = y.get_offset_element(2-1);
            confuseus[15-1] = x.get_offset_element(2-1);
            confuseus[16-1] = y.get_offset_element(2-1);
        }

        else
        {
            confuseus[8-1]  = y.get_offset_element(2-1);
            confuseus[9-1]  = x.get_offset_element(2-1);
            confuseus[15-1] = x.get_offset_element(2-1);
            confuseus[16-1] = y.get_offset_element(2-1);
        }
    }

    else
    {
        if ( x.get_effective_size() <= y.get_effective_size() )
        {
            confuseus[8-1]  = (double) x.get_effective_size();
            confuseus[9-1]  = (double) y.get_effective_size();
            confuseus[15-1] = (double) x.get_effective_size();
            confuseus[16-1] = (double) y.get_effective_size();
        }

        else
        {
            confuseus[8-1]  = (double) y.get_effective_size();
            confuseus[9-1]  = (double) x.get_effective_size();
            confuseus[15-1] = (double) x.get_effective_size();
            confuseus[16-1] = (double) y.get_effective_size();
        }
    }

    if ( !is_cross_kernel_valid )
    {
        if ( sparse_kernel )
        {
            calc_cross_kernel((long) x.get_offset_element(2-1));
        }

        else
        {
            calc_cross_kernel(x.get_effective_size());
        }
    }

    try
    {
        result = cross_kernel_def(&DoubleArgEvaluationFunctionKernel,(void *) argumentus);
    }

    catch ( int x )
    {
        L_THROW(x);
    }

    if ( !gsl_finite(result) )
    {
        L_THROW(0);
    }

    return result;
}

basic_kernel basic_kernel::kernel_deriv(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz)
{
    void *argumentus[11];
    double confuseus[16];
    basic_kernel result(*this);
    long i,j,k;

    argumentus[1-1]  = NULL;
    argumentus[2-1]  = NULL;
    argumentus[3-1]  = (void *) &uv;
    argumentus[4-1]  = (void *) &uc;
    argumentus[5-1]  = (void *) confuseus;
    argumentus[6-1]  = (void *) &x;
    argumentus[7-1]  = (void *) &y;
    argumentus[8-1]  = (void *) &mean;
    argumentus[9-1]  = NULL;
    argumentus[10-1] = (void *) &covar;
    argumentus[11-1] = (void *) &scale_prod;

    confuseus[1-1]  = (double) xd;
    confuseus[2-1]  = (double) yd;
    confuseus[3-1]  = (double) xi;
    confuseus[4-1]  = (double) yi;
    confuseus[5-1]  = (double) xz;
    confuseus[6-1]  = (double) yz;
    confuseus[7-1]  = (double) covw;

    confuseus[10-1] = GSL_NAN;
    confuseus[11-1] = GSL_NAN;
    confuseus[12-1] = GSL_NAN;
    confuseus[13-1] = GSL_NAN;
    confuseus[14-1] = GSL_NAN;

    if ( sparse_kernel )
    {
        if ( x.get_offset_element(2-1) <= y.get_offset_element(2-1) )
        {
            confuseus[8-1]  = x.get_offset_element(2-1);
            confuseus[9-1]  = y.get_offset_element(2-1);
            confuseus[15-1] = x.get_offset_element(2-1);
            confuseus[16-1] = y.get_offset_element(2-1);
        }

        else
        {
            confuseus[8-1]  = y.get_offset_element(2-1);
            confuseus[9-1]  = x.get_offset_element(2-1);
            confuseus[15-1] = x.get_offset_element(2-1);
            confuseus[16-1] = y.get_offset_element(2-1);
        }
    }

    else
    {
        if ( x.get_effective_size() <= y.get_effective_size() )
        {
            confuseus[8-1]  = (double) x.get_effective_size();
            confuseus[9-1]  = (double) y.get_effective_size();
            confuseus[15-1] = (double) x.get_effective_size();
            confuseus[16-1] = (double) y.get_effective_size();
        }

        else
        {
            confuseus[8-1]  = (double) y.get_effective_size();
            confuseus[9-1]  = (double) x.get_effective_size();
            confuseus[15-1] = (double) x.get_effective_size();
            confuseus[16-1] = (double) y.get_effective_size();
        }
    }

    if ( !is_kernel_deriv_valid )
    {
        calc_kernel_deriv();
    }

    if ( ( k = uv.get_effective_size() ) > 0 )
    {
        for ( i = 1 ; i <= k ; i++ )
        {
            confuseus[10-1] = 3;
            confuseus[11-1] = i;

            try
            {
                (result.uv)[i-1] = kernel_deriv_def(&DoubleArgEvaluationFunctionKernel,(void *) argumentus);
            }

            catch ( int x )
            {
                L_THROW(x);
            }

            if ( !gsl_finite((result.uv)[i-1]) )
            {
                L_THROW(0);
            }
        }
    }

    if ( !fix_covw )
    {
        confuseus[10-1] = 5;
        confuseus[11-1] = 7;

        try
        {
            result.covw = kernel_deriv_def(&DoubleArgEvaluationFunctionKernel,(void *) argumentus);
        }

        catch ( int x )
        {
            L_THROW(x);
        }

        if ( !gsl_finite(result.covw) )
        {
            L_THROW(0);
        }
    }

    if ( !fix_mean )
    {
        if ( ( k = mean.get_effective_size() ) > 0 )
        {
            for ( i = 1 ; i <= k ; i++ )
            {
                confuseus[10-1] = 8;
                confuseus[11-1] = i;

                try
                {
                    (result.mean)[i-1] = kernel_deriv_def(&DoubleArgEvaluationFunctionKernel,(void *) argumentus);
                }

                catch ( int x )
                {
                    L_THROW(x);
                }

                if ( !gsl_finite((result.mean)[i-1]) )
                {
                    L_THROW(0);
                }
            }
        }
    }

    if ( !fix_covar )
    {
        /* Assumption - the matrix is square */

        if ( ( k = covar.get_effective_height() ) > 0 )
        {
            for ( i = 1 ; i <= k ; i++ )
            {
                for ( j = 1 ; j <= i ; j++ )
                {
                    confuseus[10-1] = 9+i;
                    confuseus[11-1] = j;

                    try
                    {
                        (result.covar)(i-1,j-1) = kernel_deriv_def(&DoubleArgEvaluationFunctionKernel,(void *) argumentus);
                    }

                    catch ( int x )
                    {
                        L_THROW(x);
                    }

                    if ( !gsl_finite((result.covar)(i-1,j-1)) )
                    {
                        L_THROW(0);
                    }
                }
            }
        }
    }

    return result;
}




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

void basic_kernel::scale_kernel(double scale)
{
    uv *= scale;

    if ( !fix_covw )
    {
        covw *= scale;
    }

    if ( !fix_mean )
    {
        mean *= scale;
    }

    if ( !fix_covar )
    {
        covar *= scale;
    }

    return;
}

void basic_kernel::add_kernel(const basic_kernel &what, double scale)
{
    is_opt            = 0;
    optim_kernel      = NULL;
    optim_kernel_fast = NULL;
    
    uv += ( what.uv * scale );

    if ( !fix_covw )
    {
        covw += ( what.covw * scale );
    }

    if ( !fix_mean )
    {
        mean += ( what.mean * scale );
    }

    if ( !fix_covar )
    {
        covar += ( what.covar * scale );
    }

    return;
}




//
// Information
// ===========
//

int basic_kernel::is_sparse(void) const
{
    return sparse_kernel;
}




//
// IO stuff
// ========
//

void basic_kernel::prime_io(std::ostream *_where_to, int _echo_level)
{
    where_to   = _where_to;
    echo_level = _echo_level;

    return;
}




//
// Kernel control functions
// ========================
//

void basic_kernel::simplify_kernel(void)
{
    try
    {
        if ( sparse_kernel )
        {
            kernel_def.naive_replace(DOT_PRODUCT_SPARSE,5,12);
            kernel_def.naive_replace(DIFF_MAG_SQUARED_SPARSE,5,13);
            kernel_def.naive_replace(DIFF_MAGNITUDE_SPARSE,5,14);

            /*
               Mean and covariance correction is not possible for
               sparse vectors.
            */
        }

        else
        {
            kernel_def.naive_replace(DOT_PRODUCT_NONSPARSE,5,12);
            kernel_def.naive_replace(DIFF_MAG_SQUARED_NONSPARSE,5,13);
            kernel_def.naive_replace(DIFF_MAGNITUDE_NONSPARSE,5,14);

            if ( fix_covar )
            {
                if ( covar.is_diagonal() || covar.is_scalar() )
                {
                    /*
                       We can safely use the diagonal form here, as if the
                       covar matrix is fixed then the diagonals will remain
                       0 permanently.
                    */

                    kernel_def.naive_replaceRow(VAR_REP_FORM_X_DIAG,1);
                    kernel_def.naive_replaceRow(VAR_REP_FORM_Y_DIAG,2);
                }

                else
                {
                    kernel_def.naive_replaceRow(VAR_REP_FORM_X_GENERAL,1);
                    kernel_def.naive_replaceRow(VAR_REP_FORM_Y_GENERAL,2);
                }
            }

            else
            {
                kernel_def.naive_replaceRow(VAR_REP_FORM_X_GENERAL,1);
                kernel_def.naive_replaceRow(VAR_REP_FORM_Y_GENERAL,2);
            }
        }
    }

    catch ( int x )
    {
        L_THROW(x);
    }

    return;
}

void basic_kernel::calc_cross_kernel(long dim)
{
    long i;

    THROW_ASSERT(dim >= 0);

    try
    {
        if ( !is_cross_kernel_valid )
        {
            cross_kernel_def = kernel_def;

            /*
               We must assume that the dimension will not change (if it does
               then this equation will be invalid.  But then, unless dim is
               fixed, the whole process will be flawed).  Also, the
               simplification process can work more efficiently if the some
               arguments are removed.  We can remove:

               var(5,8)  = dim
               var(5,9)  = dim
               var(5,15) = dim
               var(5,16) = dim
            */

            cross_kernel_def.fast_naive_replace(dim,5,8);
            cross_kernel_def.fast_naive_replace(dim,5,9);
            cross_kernel_def.fast_naive_replace(dim,5,15);
            cross_kernel_def.fast_naive_replace(dim,5,16);

            cross_kernel_def.simplify();

            if ( dim > 0 )
            {
                for ( i = 1 ; i <= dim ; i++ )
                {
                    cross_kernel_def = cross_kernel_def.diff(6,i);
                }
            }

            is_cross_kernel_valid = 1;
        }
    }

    catch ( int x )
    {
        L_THROW(x);
    }

    return;
}

void basic_kernel::calc_kernel_deriv(void)
{
    try
    {
        if ( !is_kernel_deriv_valid )
        {
            is_kernel_deriv_valid = 1;

            kernel_deriv_def = kernel_def;

            kernel_deriv_def = kernel_deriv_def.diffMatrix(5,10,5,11);
        }
    }

    catch ( int x )
    {
        L_THROW(x);
    }

    return;
}




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

double DoubleArgEvaluationFunctionKernel(void *XargContents, long i, long j)
{
    double result;
    void **argContents;

    argContents = (void **) XargContents;

    switch ( i )
    {
        case 3:
        {
            result = (double) ((*((fVECTOR *) argContents[3-1]))[j-1]);

            break;
        }

        case 4:
        {
            result = (double) ((*((fVECTOR *) argContents[4-1]))[j-1]);

            break;
        }

        case 5:
        {
            if ( ( j >= 1 ) && ( j <= 16 ) )
            {
                result = ((double *) argContents[5-1])[j-1];
            }

            else
            {
                result = GSL_NAN;
            }

            break;
        }

        case 6:
        {
            result = (double) ((*((fVECTOR *) argContents[6-1]))[j-1]);

            if ( ((int *) argContents[11-1])[0] )
            {
                result /= ((double *) argContents[5-1])[8-1];
            }

            break;
        }

        case 7:
        {
            result = (double) ((*((fVECTOR *) argContents[7-1]))[j-1]);

            if ( ((int *) argContents[11-1])[0] )
            {
                result /= ((double *) argContents[5-1])[8-1];
            }

            break;
        }

        case 8:
        {
            result = (double) ((*((fVECTOR *) argContents[8-1]))[j-1]);

            break;
        }

        default:
        {
            if ( i >= 10 )
            {
                result = (double) ((*((fMATRIX *) argContents[10-1]))(i-10,j-1));
            }

            else
            {
                result = GSL_NAN;
            }

            break;
        }
    }

    return result;
}




//
// Stream io
// =========
//

std::ostream &operator<<(std::ostream &output, const basic_kernel &dumpee)
{
    try
    {
        output << "Kernel: "       << dumpee.kernel_def       << "\n";
        output << "Cross kernel: " << dumpee.cross_kernel_def << "\n";
        output << "Kernel deriv: " << dumpee.kernel_deriv_def << "\n\n";

        output << "Cross kernel validity: " << dumpee.is_cross_kernel_valid << "\n";
        output << "Kernel deriv validity: " << dumpee.is_kernel_deriv_valid << "\n\n";

        output << "Scale products: " << dumpee.scale_prod    << "\n";
        output << "Sparse kernel:  " << dumpee.sparse_kernel << "\n\n";

        output << "uc: " << dumpee.uc << "\n";
        output << "uv: " << dumpee.uv << "\n\n";

        output << "covw: "  << fixnum(dumpee.covw) << "\n";
        output << "mean: "  << dumpee.mean         << "\n";
        output << "covar: " << dumpee.covar        << "\n\n";

        output << "fix_mean: "  << dumpee.fix_mean  << "\n";
        output << "fix_covar: " << dumpee.fix_covar << "\n";
        output << "fix_covw: "  << dumpee.fix_covw  << "\n\n";

        output << "is_opt: "  << dumpee.is_opt  << "\n\n";
    }

    catch ( int x )
    {
        L_THROW(x);
    }

    return output;
}

std::istream &operator>>(std::istream &input, basic_kernel &dest)
{
    long i,kerauto;
    char selauto;
    wait_dummy howzat;

    try
    {
        if ( dest.where_to == NULL )
        {
            input >> howzat; input >> dest.kernel_def;
            input >> howzat; input >> dest.cross_kernel_def;
            input >> howzat; input >> dest.kernel_deriv_def;

            input >> howzat; input >> dest.is_cross_kernel_valid;
            input >> howzat; input >> dest.is_kernel_deriv_valid;

            input >> howzat; input >> dest.scale_prod;
            input >> howzat; input >> dest.sparse_kernel;

            input >> howzat; input >> dest.uc;
            input >> howzat; input >> dest.uv;

            input >> howzat; input >> dest.covw;
            input >> howzat; input >> dest.mean;
            input >> howzat; input >> dest.covar;

            input >> howzat; input >> dest.fix_mean;
            input >> howzat; input >> dest.fix_covar;
            input >> howzat; input >> dest.fix_covw;

            input >> howzat; input >> dest.is_opt;

            kerauto = dest.is_opt;
        }

        else
        {
            dest.is_opt = 0;
            kerauto = 0;

            dest.cross_kernel_def = "0";
            dest.kernel_deriv_def = "0";

            dest.is_cross_kernel_valid = 0;
            dest.is_kernel_deriv_valid = 0;

            *(dest.where_to) << "Use inbuilt kernel (y = yes): ";
            input >> selauto;
            if ( dest.echo_level ) { *(dest.where_to) << selauto << "\n"; }

            if ( selauto == 'y' )
            {
                repeater:

                *(dest.where_to) << "Kernel number (0 for details): ";
                input >> kerauto;
                if ( dest.echo_level ) { *(dest.where_to) << kerauto << "\n"; }

                if ( ( kerauto >= 1 ) && ( kerauto <= N_KERNS ) )
                {
                    *(dest.where_to) << "Kernel: ";
                    *(dest.where_to) << inbuilt_kernels[kerauto-1] << "\n";

                    dest.kernel_def = inbuilt_kernels[kerauto-1];
                }

                else
                {
                    for ( kerauto = 1 ; kerauto <= N_KERNS ; kerauto++ )
                    {
                        *(dest.where_to) << kern_descr[kerauto-1];
                    }

                    *(dest.where_to) << "Where: - non-sparse indicates that this kernel will not work with\n";
                    *(dest.where_to) << "         sparse vectors.\n";
                    *(dest.where_to) << "       - inbuilt indicates that the kernel function has been hardcoded\n";
                    *(dest.where_to) << "         for fast excecution (non inbuilt kernels may be prohibitively\n";
                    *(dest.where_to) << "         slow if frequent kernel operations are required - eg. SMO on\n";
                    *(dest.where_to) << "         large dataset).\n";
                    *(dest.where_to) << "       - z_i = x_i - y_i\n";

                    goto repeater;
                }
            }

            else
            {
                *(dest.where_to) << "Kernel: ";
                input >> dest.kernel_def;
                if ( dest.echo_level ) { *(dest.where_to) << dest.kernel_def << "\n"; }
            }

            if ( selauto == 'y' )
            {
                if ( inbuilt_var_use[kerauto-1][0] > 0 )
                {
                    (dest.uc).make_zero(DO_FORCE);
                    (dest.uc).make_normal(DO_FORCE);
                    (dest.uc).pad_vector(inbuilt_var_use[kerauto-1][0]);

                    *(dest.where_to) << "User constants (" << inbuilt_var_use[kerauto-1][0] << "):\n";

                    for ( i = 1 ; i <= inbuilt_var_use[kerauto-1][0] ; i++ )
                    {
                        *(dest.where_to) << "var(4," << i << ") = ";
                        input >> (dest.uc)[i-1];
                        if ( dest.echo_level ) { *(dest.where_to) << (dest.uc)[i-1] << "\n"; }

                        (dest.kernel_def).naive_replace((dest.uc)[i-1],4,i);
                    }
                }

                if ( inbuilt_var_use[kerauto-1][1] > 0 )
                {
                    (dest.uv).make_zero(DO_FORCE);
                    (dest.uv).make_normal(DO_FORCE);
                    (dest.uv).pad_vector(inbuilt_var_use[kerauto-1][1]);

                    *(dest.where_to) << "User variables (" << inbuilt_var_use[kerauto-1][1] << "):\n";

                    for ( i = 1 ; i <= inbuilt_var_use[kerauto-1][1] ; i++ )
                    {
                        *(dest.where_to) << "var(3," << i << ") = ";
                        input >> (dest.uv)[i-1];
                        if ( dest.echo_level ) { *(dest.where_to) << (dest.uv)[i-1] << "\n"; }
                    }
                }
            }

            else
            {
                *(dest.where_to) << "User constants:\n";
                (dest.uc).prime_io(dest.where_to,dest.echo_level);
                input >> dest.uc;

                *(dest.where_to) << "User variables:\n";
                (dest.uv).prime_io(dest.where_to,dest.echo_level);
                input >> dest.uv;
            }

            *(dest.where_to) << "Scale products (0 = no): ";
            input >> dest.scale_prod;
            if ( dest.echo_level ) { *(dest.where_to) << dest.scale_prod << "\n"; }

            *(dest.where_to) << "Sparse kernel (0 = no, 1 = yes ): ";
            input >> dest.sparse_kernel;
            if ( dest.echo_level ) { *(dest.where_to) << dest.sparse_kernel << "\n"; }

            if ( dest.sparse_kernel )
            {
                (dest.mean).make_zero(DO_FORCE);
                (dest.covar).make_zero(DO_FORCE);

                dest.fix_mean  = 1;
                dest.fix_covar = 1;

                *(dest.where_to) << "Scalar variance: ";
                input >> dest.covw;
                if ( dest.echo_level ) { *(dest.where_to) << dest.covw << "\n"; }

                *(dest.where_to) << "Fix scalar variance (0 = no): ";
                input >> dest.fix_covw;
                if ( dest.echo_level ) { *(dest.where_to) << dest.fix_covw << "\n"; }
            }

            else
            {
                *(dest.where_to) << "Scalar variance: ";
                input >> dest.covw;
                if ( dest.echo_level ) { *(dest.where_to) << dest.covw << "\n"; }

                *(dest.where_to) << "Mean correction:\n";
                (dest.mean).prime_io(dest.where_to,dest.echo_level);
                input >> dest.mean;

                *(dest.where_to) << "Matrix covariance correction (lower triangular):\n";
                (dest.covar).prime_io(dest.where_to,dest.echo_level);
                input >> dest.covar;

                *(dest.where_to) << "Fix scalar variance (0 = no): ";
                input >> dest.fix_covw;
                if ( dest.echo_level ) { *(dest.where_to) << dest.fix_covw << "\n"; }

                *(dest.where_to) << "Fix mean correction (0 = no): ";
                input >> dest.fix_mean;
                if ( dest.echo_level ) { *(dest.where_to) << dest.fix_mean << "\n"; }

                *(dest.where_to) << "Fix matrix covariance correction (0 = no): ";
                input >> dest.fix_covar;
                if ( dest.echo_level ) { *(dest.where_to) << dest.fix_covar << "\n"; }
            }

            if ( kerauto )
            {
                if ( dest.sparse_kernel )
                {
                    if ( (inbuilt_opt_kernel[kerauto-1]).kernel_sparse != NULL )
                    {
                        dest.is_opt = kerauto;
                    }
                }

                else
                {
                    if ( dest.scale_prod )
                    {
                        if ( (inbuilt_opt_kernel[kerauto-1]).kernel_scale_prod != NULL )
                        {
                            dest.is_opt = kerauto;
                        }
                    }

                    else
                    {
                        if ( (inbuilt_opt_kernel[kerauto-1]).kernel != NULL )
                        {
                            dest.is_opt = kerauto;
                        }
                    }
                }
            }

            dest.simplify_kernel();
            dest.prime_io(NULL,0);
        }

        dest.optim_kernel      = NULL;
        dest.optim_kernel_fast = NULL;

        if ( dest.is_opt )
        {
            if ( dest.sparse_kernel )
            {
                dest.optim_kernel      = (inbuilt_opt_kernel[kerauto-1]).kernel_sparse;
                dest.optim_kernel_fast = (inbuilt_opt_kernel[kerauto-1]).kernel_sparse_opt;
            }

            else
            {
                if ( dest.scale_prod )
                {
                    dest.optim_kernel      = (inbuilt_opt_kernel[kerauto-1]).kernel_scale_prod;
                    dest.optim_kernel_fast = NULL;
                }

                else
                {
                    dest.optim_kernel      = (inbuilt_opt_kernel[kerauto-1]).kernel;
                    dest.optim_kernel_fast = NULL;
                }
            }
        }
    }

    catch ( int x )
    {
        L_THROW(x);
    }

    return input;
}




//
// Operator overloading
// ====================
//
// These are just like the kernel manipulation functions described
// previously - the same provisos apply.
//
// +  posation                  - unary,  return rvalue
// +  addition                  - binary, return rvalue
// += additive assignment       - binary, return lvalue
// *  multiplication            - binary, return rvalue
// /  division                  - binary, return rvalue
// *= multiplicative assignment - binary, return lvalue
// /= divisive       assignment - binary, return lvalue
//

basic_kernel operator+(const basic_kernel &left_op)
{
    basic_kernel result(left_op);

    return result;
}

basic_kernel operator-(const basic_kernel &left_op)
{
    basic_kernel result(left_op);

    result.scale_kernel(-1.0);

    return result;
}

basic_kernel operator+(const basic_kernel &left_op, const basic_kernel &right_op)
{
    basic_kernel result(left_op);

    result.add_kernel(right_op,1.0);

    return result;
}

basic_kernel operator-(const basic_kernel &left_op, const basic_kernel &right_op)
{
    basic_kernel result(left_op);

    result.add_kernel(right_op,-1.0);

    return result;
}

basic_kernel operator*(const basic_kernel &left_op, const double &right_op)
{
    basic_kernel result(left_op);

    result.scale_kernel(right_op);

    return result;
}

basic_kernel operator*(const double &right_op, const basic_kernel &left_op)
{
    basic_kernel result(left_op);

    result.scale_kernel(right_op);

    return result;
}

basic_kernel operator/(const basic_kernel &left_op, const double &right_op)
{
    basic_kernel result(left_op);

    result.scale_kernel(1/right_op);

    return result;
}

basic_kernel &operator+=(basic_kernel &left_op, const basic_kernel &right_op)
{
    left_op.add_kernel(right_op,1.0);

    return left_op;
}

basic_kernel &operator-=(basic_kernel &left_op, const basic_kernel &right_op)
{
    left_op.add_kernel(right_op,-1.0);

    return left_op;
}

basic_kernel &operator*=(basic_kernel &left_op, const double &right_op)
{
    left_op.scale_kernel(right_op);

    return left_op;
}

basic_kernel &operator/=(basic_kernel &left_op, const double &right_op)
{
    left_op.scale_kernel(1.0/right_op);

    return left_op;
}




/******************************************************************/
/*                                                                */
/* kernel dictionary class.                                       */
/*                                                                */
/******************************************************************/

kernel_dicto::kernel_dicto()
{
    where_to = NULL;
    echo_level = 0;

    dsize       = 1;
    fix_weights = 1;

    REPNEWB(kern,basic_kernel *,dsize);

    kweights.make_zero(DO_FORCE);
    kweights.make_normal(DO_FORCE);

    REPNEW((kern[0]),basic_kernel);
    kweights.addend(1.0);

    return;
}

kernel_dicto::kernel_dicto(const kernel_dicto &what)
{
    long i;

    where_to = NULL;
    echo_level = 0;

    dsize       = what.dsize;
    fix_weights = what.fix_weights;

    THROW_ASSERT(dsize > 0);

    REPNEWB(kern,basic_kernel *,dsize);

    kweights.make_zero(DO_FORCE);
    kweights.make_normal(DO_FORCE);

    for ( i = 0 ; i < dsize ; i++ )
    {
        REPNEW((kern[i]),basic_kernel((*(what.kern)[i])));
        kweights.addend((what.kweights).get_offset_element(i));
    }

    return;
}

kernel_dicto::kernel_dicto(long _dsize, int _fix_weights, double *weights, char **kern_def, fVECTOR **_uc, fVECTOR **_uv, double *_covw, fVECTOR **_mean, fMATRIX **_covar, int *_fix_covw, int *_fix_mean, int *_fix_covar, int *_scale_prod, int _sparse_kernel)
{
    long i;

    THROW_ASSERT(_dsize > 0);

    where_to = NULL;
    echo_level = 0;

    dsize       = _dsize;
    fix_weights = _fix_weights;

    REPNEWB(kern,basic_kernel *,dsize);

    kweights.make_zero(DO_FORCE);
    kweights.make_normal(DO_FORCE);

    for ( i = 0 ; i < dsize ; i++ )
    {
        REPNEW((kern[i]),basic_kernel(kern_def[i],*(_uc[i]),*(_uv[i]),_covw[i],*(_mean[i]),*(_covar[i]),_fix_covw[i],_fix_mean[i],_fix_covar[i],_scale_prod[i],_sparse_kernel));
        kweights.addend(weights[i]);
    }

    return;
}

kernel_dicto::kernel_dicto(long _dsize, int _fix_weights, double *weights, int *kern_number, fVECTOR **_uc, fVECTOR **_uv, double *_covw, fVECTOR **_mean, fMATRIX **_covar, int *_fix_covw, int *_fix_mean, int *_fix_covar, int *_scale_prod, int _sparse_kernel)
{
    long i;

    THROW_ASSERT(_dsize > 0);

    where_to = NULL;
    echo_level = 0;

    dsize       = _dsize;
    fix_weights = _fix_weights;

    REPNEWB(kern,basic_kernel *,dsize);

    kweights.make_zero(DO_FORCE);
    kweights.make_normal(DO_FORCE);

    for ( i = 0 ; i < dsize ; i++ )
    {
        REPNEW((kern[i]),basic_kernel(kern_number[i],*(_uc[i]),*(_uv[i]),_covw[i],*(_mean[i]),*(_covar[i]),_fix_covw[i],_fix_mean[i],_fix_covar[i],_scale_prod[i],_sparse_kernel));
        kweights.addend(weights[i]);
    }

    return;
}

kernel_dicto::kernel_dicto(double weight, int kern_number, const fVECTOR &_uc, const fVECTOR &_uv, double _covw, const fVECTOR &_mean, const fMATRIX &_covar, int _fix_covw, int _fix_mean, int _fix_covar, int _scale_prod, int _sparse_kernel)
{
    where_to = NULL;
    echo_level = 0;

    dsize       = 1;
    fix_weights = 1;

    REPNEWB(kern,basic_kernel *,dsize);

    kweights.make_zero(DO_FORCE);
    kweights.make_normal(DO_FORCE);

    REPNEW((kern[0]),basic_kernel(kern_number,_uc,_uv,_covw,_mean,_covar,_fix_covw,_fix_mean,_fix_covar,_scale_prod,_sparse_kernel));
    kweights.addend(weight);

    return;
}

kernel_dicto &kernel_dicto::operator=(const kernel_dicto &source)
{
    long i;

    if ( this != &source )
    {
        for ( i = 0 ; i < dsize ; i++ )
        {
            REPDEL((kern[i]));
        }

        REPDELB(kern);

        dsize       = source.dsize;
        fix_weights = source.fix_weights;

        THROW_ASSERT(dsize > 0);

        REPNEWB(kern,basic_kernel *,dsize);

        kweights.make_zero(DO_FORCE);
        kweights.make_normal(DO_FORCE);

        for ( i = 0 ; i < dsize ; i++ )
        {
            REPNEW((kern[i]),basic_kernel((*(source.kern)[i])));
            kweights.addend((source.kweights).get_offset_element(i));
        }
    }

    return *this;
}

kernel_dicto::~kernel_dicto()
{
    long i;

    THROW_ASSERT(dsize > 0);

    for ( i = 0 ; i < dsize ; i++ )
    {
        REPDEL((kern[i]));
    }

    REPDELB(kern);

    return;
}




//
// Normalisation
// =============
//
// make kweights'.kweights = normal
//

void kernel_dicto::normalise(double normal)
{
    double totweight;

    totweight = sqrt((double) (kweights*kweights));

    kweights /= totweight;
    kweights *= normal;

    return;
}




//
// Kernel evaluation functions
// ===========================
//
// Calculate the relevant result.
//

double kernel_dicto::kernel(const fVECTOR &x, const fVECTOR &y, L_DOUBLE *x_fast, L_DOUBLE *y_fast, double xd, double yd, long xi, long yi, double xz, double yz)
{
    #ifdef DO__FLOPS
    REG_kernel_calls;
    #endif

    double result = 0;
    long i;

    for ( i = 0 ; i < dsize ; i++ )
    {
        result += ((double) kweights[i]) * ((kern[i])->kernel(x,y,x_fast,y_fast,xd,yd,xi,yi,xz,yz));
    }

    return result;
}


double kernel_dicto::cross_kernel(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz)
{
    #ifdef DO__FLOPS
    REG_c_kernel_calls;
    #endif

    double result = 0;
    long i;

    for ( i = 0 ; i < dsize ; i++ )
    {
        result += ((double) kweights[i]) * ((kern[i])->cross_kernel(x,y,xd,yd,xi,yi,xz,yz));
    }

    return result;
}


kernel_dicto kernel_dicto::kernel_deriv(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz)
{
    #ifdef DO__FLOPS
    REG_d_kernel_calls;
    #endif

    kernel_dicto result(*this);
    long i;

    for ( i = 0 ; i < dsize ; i++ )
    {
        if ( !fix_weights )
        {
            (result.kweights)[i] = (kern[i])->kernel(x,y,NULL,NULL,xd,yd,xi,yi,xz,yz);
        }

        else
        {
            (result.kweights)[i] = 0;
        }

        (*((result.kern)[i])) = (kern[i])->kernel_deriv(x,y,xd,yd,xi,yi,xz,yz);
    }

    return result;
}




//
// Elementwise access functions
// ============================
//

basic_kernel &kernel_dicto::kernel_element(long i)
{
    THROW_ASSERT(i>=1);
    THROW_ASSERT(i<=dsize);

    return *(kern[i-1]);
}

L_DOUBLE &kernel_dicto::weight_element(long i)
{
    THROW_ASSERT(i>=1);
    THROW_ASSERT(i<=dsize);

    return kweights[i-1];
}

long kernel_dicto::get_dsize(void)
{
    return dsize;
}

fVECTOR &kernel_dicto::getweights(void)
{
    return kweights;
}




//
// Kernel manipulation functions
// =============================
//

void kernel_dicto::scale_kernel(double scale)
{
    long i;

    if ( !fix_weights )
    {
        kweights *= scale;
    }

    for ( i = 0 ; i < dsize ; i++ )
    {
        *(kern[i]) *= scale;
    }
}

void kernel_dicto::add_kernel(const kernel_dicto &what, double scale)
{
    long i;

    if ( !fix_weights )
    {
        kweights += ( (what.kweights) * scale );
    }

    for ( i = 0 ; i < dsize ; i++ )
    {
        *(kern[i]) += ( (*((what.kern)[i])) * scale );
    }
}




//
// Information
// ===========
//

int kernel_dicto::is_sparse(void) const
{
    return (kern[0])->is_sparse();
}




//
// IO stuff
// ========
//

void kernel_dicto::prime_io(std::ostream *_where_to, int _echo_level)
{
    where_to   = _where_to;
    echo_level = _echo_level;

    return;
}




//
// Stream io
// =========
//

std::ostream &operator<<(std::ostream &output, const kernel_dicto &dumpee)
{
    long i;

    output << "dsize: "       << dumpee.dsize       << "\n";
    output << "fix_weights: " << dumpee.fix_weights << "\n\n";

    for ( i = 0 ; i < dumpee.dsize ; i++ )
    {
        output << "Weight " << i+1 << ": "  << fixnum((dumpee.kweights).get_offset_element(i)) << "\n";
        output << "Kernel " << i+1 << ":\n" << *((dumpee.kern)[i])                             << "\n\n";
    }

    return output;
}

std::istream &operator>>(std::istream &input, kernel_dicto &dest)
{
    long i;
    double x;
    wait_dummy howzat;

    for ( i = 0 ; i < dest.dsize ; i++ )
    {
        REPDEL(((dest.kern)[i]));
    }

    REPDELB((dest.kern));

    if ( dest.where_to == NULL )
    {
        input >> howzat; input >> dest.dsize;
        input >> howzat; input >> dest.fix_weights;

        REPNEWB((dest.kern),basic_kernel *,(dest.dsize));

        (dest.kweights).make_zero(DO_FORCE);
        (dest.kweights).make_normal(DO_FORCE);

        for ( i = 0 ; i < dest.dsize ; i++ )
        {
            REPNEW(((dest.kern)[i]),basic_kernel);

            input >> howzat; input >> x;

            (dest.kweights).addend(x);

            input >> howzat; input >> *((dest.kern)[i]);
        }
    }

    else
    {
        *(dest.where_to) << "Dictionary size: ";
        input >> dest.dsize;
        if ( dest.echo_level ) { *(dest.where_to) << dest.dsize << "\n"; }

        REPNEWB((dest.kern),basic_kernel *,(dest.dsize));

        (dest.kweights).make_zero(DO_FORCE);
        (dest.kweights).make_normal(DO_FORCE);

        for ( i = 0 ; i < dest.dsize ; i++ )
        {
            *(dest.where_to) << "Kernel number " << i << ".\n\n";

            REPNEW(((dest.kern)[i]),basic_kernel);

            ((dest.kern)[i])->prime_io(dest.where_to,dest.echo_level);
            input >> (*((dest.kern)[i]));

            *(dest.where_to) << "Kernel weight: ";
            input >> x;
            if ( dest.echo_level ) { *(dest.where_to) << x << "\n"; }

            (dest.kweights).addend(x);

            *(dest.where_to) << "\n";
        }

        *(dest.where_to) << "Fix weights (0 = no): ";
        input >> dest.fix_weights;
        if ( dest.echo_level ) { *(dest.where_to) << dest.fix_weights << "\n"; }

        dest.prime_io(NULL,0);
    }

    return input;
}




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

kernel_dicto operator+(const kernel_dicto &left_op)
{
    kernel_dicto result(left_op);

    return result;
}

kernel_dicto operator-(const kernel_dicto &left_op)
{
    kernel_dicto result(left_op);
    long i;

    for ( i = 0 ; i < result.dsize ; i++ )
    {
        *((result.kern)[i]) *= -1.0;
    }

    result.kweights *= -1.0;

    return result;
}

kernel_dicto operator+(const kernel_dicto &left_op, const kernel_dicto &right_op)
{
    kernel_dicto result(left_op);
    long i;

    THROW_ASSERT(left_op.dsize == right_op.dsize);

    for ( i = 0 ; i < result.dsize ; i++ )
    {
        *((result.kern)[i]) += *((right_op.kern)[i]);
    }

    result.kweights += right_op.kweights;

    return result;
}

kernel_dicto operator-(const kernel_dicto &left_op, const kernel_dicto &right_op)
{
    kernel_dicto result(left_op);
    long i;

    THROW_ASSERT(left_op.dsize == right_op.dsize);

    for ( i = 0 ; i < result.dsize ; i++ )
    {
        *((result.kern)[i]) -= *((right_op.kern)[i]);
    }

    result.kweights -= right_op.kweights;

    return result;
}

kernel_dicto operator*(const kernel_dicto &left_op, const double &right_op)
{
    kernel_dicto result(left_op);
    long i;

    for ( i = 0 ; i < result.dsize ; i++ )
    {
        *((result.kern)[i]) *= right_op;
    }

    result.kweights *= right_op;

    return result;
}

kernel_dicto operator*(const double &right_op, const kernel_dicto &left_op)
{
    kernel_dicto result(left_op);
    long i;

    for ( i = 0 ; i < result.dsize ; i++ )
    {
        *((result.kern)[i]) *= right_op;
    }

    result.kweights *= right_op;

    return result;
}

kernel_dicto operator/(const kernel_dicto &left_op, const double &right_op)
{
    kernel_dicto result(left_op);
    long i;

    for ( i = 0 ; i < result.dsize ; i++ )
    {
        *((result.kern)[i]) /= right_op;
    }

    result.kweights /= right_op;

    return result;
}

kernel_dicto &operator+=(kernel_dicto &left_op, const kernel_dicto &right_op)
{
    long i;

    THROW_ASSERT(left_op.dsize == right_op.dsize);

    for ( i = 0 ; i < left_op.dsize ; i++ )
    {
        *((left_op.kern)[i]) += *((right_op.kern)[i]);
    }

    left_op.kweights += right_op.kweights;

    return left_op;
}

kernel_dicto &operator-=(kernel_dicto &left_op, const kernel_dicto &right_op)
{
    long i;

    THROW_ASSERT(left_op.dsize == right_op.dsize);

    for ( i = 0 ; i < left_op.dsize ; i++ )
    {
        *((left_op.kern)[i]) -= *((right_op.kern)[i]);
    }

    left_op.kweights -= right_op.kweights;

    return left_op;
}

kernel_dicto &operator*=(kernel_dicto &left_op, const double &right_op)
{
    long i;

    for ( i = 0 ; i < left_op.dsize ; i++ )
    {
        *((left_op.kern)[i]) *= right_op;
    }

    left_op.kweights *= right_op;

    return left_op;
}

kernel_dicto &operator/=(kernel_dicto &left_op, const double &right_op)
{
    long i;

    for ( i = 0 ; i < left_op.dsize ; i++ )
    {
        *((left_op.kern)[i]) /= right_op;
    }

    left_op.kweights /= right_op;

    return left_op;
}





////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////////////                                           ///////////////
//////////////                                           ///////////////
//////////////    Optimised kernels start here.          ///////////////
//////////////                                           ///////////////
//////////////                                           ///////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//
// dot_prod: xm'ym
// dif_meas: ||xm-ym||^2 = (xm-ym)'(xm-ym)
// eucldist: ||xm-ym||
//
// Non-sparse: xm = covw.covar.(x-mean)
//             ym = covw.covar.(y-mean)
//
// Sparse: xm = x
//         ym = x
//

double dot_prod(const fVECTOR &x, const fVECTOR &y, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    STOP_FLOP_COUNT;

    double result;

    fVECTOR xm(x);
    fVECTOR ym(y);

    xm -= mean;
    ym -= mean;

    xm *= covar;
    ym *= covar;

    xm *= (L_DOUBLE) covw;
    ym *= (L_DOUBLE) covw;

    result = (double) (xm*ym);

    RESTART_FLOP_COUNT

    return result;
}

double dif_meas(const fVECTOR &x, const fVECTOR &y, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    STOP_FLOP_COUNT;

    double result;

    fVECTOR z(x);

    z -= y;
    z *= covar;
    z *= (L_DOUBLE) covw;

    result = (double) (z*z);

    RESTART_FLOP_COUNT

    return result;

    // Warning killer

    result = mean.get_offset_element(0);
}

double eucldist(const fVECTOR &x, const fVECTOR &y, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return sqrt(dif_meas(x,y,covw,mean,covar));
}

double dot_prod_scale_prod(const fVECTOR &x, const fVECTOR &y, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    STOP_FLOP_COUNT;

    double result;
    L_DOUBLE dim;

    dim = x.get_effective_size();

    fVECTOR xm(x);
    fVECTOR ym(y);

    xm -= mean;
    ym -= mean;

    xm /= dim;
    ym /= dim;

    xm *= covar;
    ym *= covar;

    xm *= (L_DOUBLE) covw;
    ym *= (L_DOUBLE) covw;

    result = (double) (xm*ym);

    RESTART_FLOP_COUNT

    return result;
}

double dif_meas_scale_prod(const fVECTOR &x, const fVECTOR &y, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    STOP_FLOP_COUNT;

    double result;
    L_DOUBLE dim;

    dim = x.get_effective_size();

    fVECTOR z(x);

    z -= y;
    z /= dim;
    z *= covar;
    z *= (L_DOUBLE) covw;

    result = (double) (z*z);

    RESTART_FLOP_COUNT

    return result;

    // Warning killer

    result = mean.get_offset_element(0);
}

double eucldist_scale_prod(const fVECTOR &x, const fVECTOR &y, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return sqrt(dif_meas_scale_prod(x,y,covw,mean,covar));
}



#define KILL_WARNINGS                                                   \
{                                                                       \
    double result;                                                      \
    long __i;                                                           \
    result = xd;                                                        \
    result = yd;                                                        \
    __i = xi;                                                           \
    __i = yi;                                                           \
    result = xz;                                                        \
    result = yz;                                                        \
    result = uc.get_offset_element(0);                                  \
    result = uv.get_offset_element(0);                                  \
    result = mean.get_offset_element(0);                                \
    result = covar.get_offset_element(0,0);                             \
}


double fast_kernel_01_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return dot_prod(x,y,covw,mean,covar);

    KILL_WARNINGS;
}

double fast_kernel_02_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return dif_meas(x,y,covw,mean,covar);

    KILL_WARNINGS;
}

double fast_kernel_03_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return eucldist(x,y,covw,mean,covar);

    KILL_WARNINGS;
}

double fast_kernel_04_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow(dot_prod(x,y,covw,mean,covar),uc.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_05_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow(dot_prod(x,y,covw,mean,covar),uv.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_06_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod(x,y,covw,mean,covar),uc.get_offset_element(0))+(uc.get_offset_element(1))),uc.get_offset_element(2));

    KILL_WARNINGS;
}

double fast_kernel_07_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod(x,y,covw,mean,covar),uv.get_offset_element(0))+(uc.get_offset_element(0))),uc.get_offset_element(1));

    KILL_WARNINGS;
}

double fast_kernel_08_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod(x,y,covw,mean,covar),uc.get_offset_element(0))+(uv.get_offset_element(0))),uc.get_offset_element(1));

    KILL_WARNINGS;
}

double fast_kernel_09_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod(x,y,covw,mean,covar),uc.get_offset_element(0))+(uc.get_offset_element(1))),uv.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_10_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod(x,y,covw,mean,covar),uc.get_offset_element(0))+(uv.get_offset_element(0))),uv.get_offset_element(1));

    KILL_WARNINGS;
}

double fast_kernel_11_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod(x,y,covw,mean,covar),uv.get_offset_element(0))+(uc.get_offset_element(0))),uv.get_offset_element(1));

    KILL_WARNINGS;
}

double fast_kernel_12_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod(x,y,covw,mean,covar),uv.get_offset_element(0))+(uv.get_offset_element(1))),uc.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_13_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod(x,y,covw,mean,covar),uv.get_offset_element(0))+(uv.get_offset_element(1))),uv.get_offset_element(2));

    KILL_WARNINGS;
}

double fast_kernel_14_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return tanh(dot_prod(x,y,covw,mean,covar));

    KILL_WARNINGS;
}

double fast_kernel_15_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return tanh(dot_prod(x,y,covw,mean,covar)+uc.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_16_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return tanh(dot_prod(x,y,covw,mean,covar)+uv.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_17_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(dot_prod(x,y,covw,mean,covar));

    KILL_WARNINGS;
}

double fast_kernel_18_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(dot_prod(x,y,covw,mean,covar))-uc.get_offset_element(0);

    KILL_WARNINGS;
}

double fast_kernel_19_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(dot_prod(x,y,covw,mean,covar))-uv.get_offset_element(0);

    KILL_WARNINGS;
}   

double fast_kernel_20_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-dif_meas(x,y,covw,mean,covar));

    KILL_WARNINGS;
}

double fast_kernel_21_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-dif_meas(x,y,covw,mean,covar))-uc.get_offset_element(0);

    KILL_WARNINGS;
}

double fast_kernel_22_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-dif_meas(x,y,covw,mean,covar))-uv.get_offset_element(0);

    KILL_WARNINGS;
}

double fast_kernel_23_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-eucldist(x,y,covw,mean,covar));

    KILL_WARNINGS;
}

double fast_kernel_24_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-eucldist(x,y,covw,mean,covar))-uc.get_offset_element(0);

    KILL_WARNINGS;
}

double fast_kernel_25_(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-eucldist(x,y,covw,mean,covar))-uv.get_offset_element(0);

    KILL_WARNINGS;
}



double fast_kernel_01_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return dot_prod_scale_prod(x,y,covw,mean,covar);

    KILL_WARNINGS;
}

double fast_kernel_02_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return dif_meas_scale_prod(x,y,covw,mean,covar);

    KILL_WARNINGS;
}

double fast_kernel_03_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return eucldist_scale_prod(x,y,covw,mean,covar);

    KILL_WARNINGS;
}

double fast_kernel_04_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow(dot_prod_scale_prod(x,y,covw,mean,covar),uc.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_05_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow(dot_prod_scale_prod(x,y,covw,mean,covar),uv.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_06_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod_scale_prod(x,y,covw,mean,covar),uc.get_offset_element(0))+(uc.get_offset_element(1))),uc.get_offset_element(2));

    KILL_WARNINGS;
}

double fast_kernel_07_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod_scale_prod(x,y,covw,mean,covar),uv.get_offset_element(0))+(uc.get_offset_element(0))),uc.get_offset_element(1));

    KILL_WARNINGS;
}

double fast_kernel_08_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod_scale_prod(x,y,covw,mean,covar),uc.get_offset_element(0))+(uv.get_offset_element(0))),uc.get_offset_element(1));

    KILL_WARNINGS;
}

double fast_kernel_09_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod_scale_prod(x,y,covw,mean,covar),uc.get_offset_element(0))+(uc.get_offset_element(1))),uv.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_10_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod_scale_prod(x,y,covw,mean,covar),uc.get_offset_element(0))+(uv.get_offset_element(0))),uv.get_offset_element(1));

    KILL_WARNINGS;
}

double fast_kernel_11_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod_scale_prod(x,y,covw,mean,covar),uv.get_offset_element(0))+(uc.get_offset_element(0))),uv.get_offset_element(1));

    KILL_WARNINGS;
}

double fast_kernel_12_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod_scale_prod(x,y,covw,mean,covar),uv.get_offset_element(0))+(uv.get_offset_element(1))),uc.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_13_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(dot_prod_scale_prod(x,y,covw,mean,covar),uv.get_offset_element(0))+(uv.get_offset_element(1))),uv.get_offset_element(2));

    KILL_WARNINGS;
}

double fast_kernel_14_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return tanh(dot_prod_scale_prod(x,y,covw,mean,covar));

    KILL_WARNINGS;
}

double fast_kernel_15_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return tanh(dot_prod_scale_prod(x,y,covw,mean,covar)+uc.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_16_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return tanh(dot_prod_scale_prod(x,y,covw,mean,covar)+uv.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_17_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(dot_prod_scale_prod(x,y,covw,mean,covar));

    KILL_WARNINGS;
}

double fast_kernel_18_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(dot_prod_scale_prod(x,y,covw,mean,covar))-uc.get_offset_element(0);

    KILL_WARNINGS;
}

double fast_kernel_19_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(dot_prod_scale_prod(x,y,covw,mean,covar))-uv.get_offset_element(0);

    KILL_WARNINGS;
}

double fast_kernel_20_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-dif_meas_scale_prod(x,y,covw,mean,covar));

    KILL_WARNINGS;
}

double fast_kernel_21_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-dif_meas_scale_prod(x,y,covw,mean,covar))-uc.get_offset_element(0);

    KILL_WARNINGS;
}

double fast_kernel_22_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-dif_meas_scale_prod(x,y,covw,mean,covar))-uv.get_offset_element(0);

    KILL_WARNINGS;
}

double fast_kernel_23_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-eucldist_scale_prod(x,y,covw,mean,covar));

    KILL_WARNINGS;
}

double fast_kernel_24_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-eucldist_scale_prod(x,y,covw,mean,covar))-uc.get_offset_element(0);

    KILL_WARNINGS;
}

double fast_kernel_25_scale_prod(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-eucldist_scale_prod(x,y,covw,mean,covar))-uv.get_offset_element(0);

    KILL_WARNINGS;
}




double fast_kernel_01_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    #ifdef KERNEL_COMP
    std::cerr << "fark me 2\n";
    #endif

    return covw*covw*dot_prod_sparse(x,y);

    KILL_WARNINGS;
}

double fast_kernel_02_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return (covw*covw*dif_meas_sparse(x,y));

    KILL_WARNINGS;
}

double fast_kernel_03_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return (covw*eucldist_sparse(x,y));

    KILL_WARNINGS;
}

double fast_kernel_04_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((covw*covw*dot_prod_sparse(x,y)),uc.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_05_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow(covw*covw*dot_prod_sparse(x,y),uv.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_06_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(covw*covw*dot_prod_sparse(x,y),uc.get_offset_element(0))+(uc.get_offset_element(1))),uc.get_offset_element(2));

    KILL_WARNINGS;
}

double fast_kernel_07_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(covw*covw*dot_prod_sparse(x,y),uv.get_offset_element(0))+(uc.get_offset_element(0))),uc.get_offset_element(1));

    KILL_WARNINGS;
}

double fast_kernel_08_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(covw*covw*dot_prod_sparse(x,y),uc.get_offset_element(0))+(uv.get_offset_element(0))),uc.get_offset_element(1));

    KILL_WARNINGS;
}

double fast_kernel_09_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(covw*covw*dot_prod_sparse(x,y),uc.get_offset_element(0))+(uc.get_offset_element(1))),uv.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_10_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(covw*covw*dot_prod_sparse(x,y),uc.get_offset_element(0))+(uv.get_offset_element(0))),uv.get_offset_element(1));

    KILL_WARNINGS;
}

double fast_kernel_11_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(covw*covw*dot_prod_sparse(x,y),uv.get_offset_element(0))+(uc.get_offset_element(0))),uv.get_offset_element(1));

    KILL_WARNINGS;
}

double fast_kernel_12_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(covw*covw*dot_prod_sparse(x,y),uv.get_offset_element(0))+(uv.get_offset_element(1))),uc.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_13_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return pow((pow(covw*covw*dot_prod_sparse(x,y),uv.get_offset_element(0))+(uv.get_offset_element(1))),uv.get_offset_element(2));

    KILL_WARNINGS;
}

double fast_kernel_14_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return tanh(covw*covw*dot_prod_sparse(x,y));

    KILL_WARNINGS;
}   

double fast_kernel_15_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return tanh((covw*covw*dot_prod_sparse(x,y))+uc.get_offset_element(0));

    KILL_WARNINGS;
}

double fast_kernel_16_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return tanh((covw*covw*dot_prod_sparse(x,y))+uv.get_offset_element(0));

    KILL_WARNINGS;
}   

double fast_kernel_17_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(covw*covw*dot_prod_sparse(x,y));

    KILL_WARNINGS;
}

double fast_kernel_18_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(covw*covw*dot_prod_sparse(x,y))-uc.get_offset_element(0);

    KILL_WARNINGS;
}

double fast_kernel_19_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(covw*covw*dot_prod_sparse(x,y))-uv.get_offset_element(0);

    KILL_WARNINGS;
}

double fast_kernel_20_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-(covw*covw*dif_meas_sparse(x,y)));

    KILL_WARNINGS;
}

double fast_kernel_21_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-(covw*covw*dif_meas_sparse(x,y)))-uc.get_offset_element(0);

    KILL_WARNINGS;
}

double fast_kernel_22_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-(covw*covw*dif_meas_sparse(x,y)))-uv.get_offset_element(0);

    KILL_WARNINGS;
}

double fast_kernel_23_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-(covw*eucldist_sparse(x,y)));

    KILL_WARNINGS;
}

double fast_kernel_24_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-(covw*eucldist_sparse(x,y)))-uc.get_offset_element(0);

    KILL_WARNINGS;
}

double fast_kernel_25_sparse(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar)
{
    return exp(-(covw*eucldist_sparse(x,y)))-uv.get_offset_element(0);

    KILL_WARNINGS;
}













#define KILL_WARNINGS_B                                                 \
{                                                                       \
    uc = NULL;                                                          \
    uv = NULL;                                                          \
}


double fast_opt_kernel_01_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    #ifdef KERNEL_COMP
    std::cerr << "fark me 1\n";
    #endif

    return (*covw)*(*covw)*dot_prod_sparse_opt(x,y);

    KILL_WARNINGS_B;
}

double fast_opt_kernel_02_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return (*covw)*(*covw)*dif_meas_sparse_opt(x,y);

    KILL_WARNINGS_B;
}

double fast_opt_kernel_03_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return (*covw)*eucldist_sparse_opt(x,y);

    KILL_WARNINGS_B;
}

double fast_opt_kernel_04_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return pow(((*covw)*(*covw)*dot_prod_sparse_opt(x,y)),uc[0]);

    KILL_WARNINGS_B;
}

double fast_opt_kernel_05_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return pow((*covw)*(*covw)*dot_prod_sparse_opt(x,y),uv[0]);

    KILL_WARNINGS_B;
}

double fast_opt_kernel_06_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return pow((pow((*covw)*(*covw)*dot_prod_sparse_opt(x,y),uc[0])+(uc[1])),uc[2]);

    KILL_WARNINGS_B;
}

double fast_opt_kernel_07_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return pow((pow((*covw)*(*covw)*dot_prod_sparse_opt(x,y),uv[0])+(uc[0])),uc[1]);

    KILL_WARNINGS_B;
}

double fast_opt_kernel_08_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return pow((pow((*covw)*(*covw)*dot_prod_sparse_opt(x,y),uc[0])+(uv[0])),uc[1]);

    KILL_WARNINGS_B;
}

double fast_opt_kernel_09_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return pow((pow((*covw)*(*covw)*dot_prod_sparse_opt(x,y),uc[0])+(uc[1])),uv[0]);

    KILL_WARNINGS_B;
}

double fast_opt_kernel_10_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return pow((pow((*covw)*(*covw)*dot_prod_sparse_opt(x,y),uc[0])+(uv[0])),uv[1]);

    KILL_WARNINGS_B;
}

double fast_opt_kernel_11_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return pow((pow((*covw)*(*covw)*dot_prod_sparse_opt(x,y),uv[0])+(uc[0])),uv[1]);

    KILL_WARNINGS_B;
}

double fast_opt_kernel_12_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return pow((pow((*covw)*(*covw)*dot_prod_sparse_opt(x,y),uv[0])+(uv[1])),uc[0]);

    KILL_WARNINGS_B;
}

double fast_opt_kernel_13_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return pow((pow((*covw)*(*covw)*dot_prod_sparse_opt(x,y),uv[0])+(uv[1])),uv[2]);

    KILL_WARNINGS_B;
}

double fast_opt_kernel_14_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return tanh((*covw)*(*covw)*dot_prod_sparse_opt(x,y));

    KILL_WARNINGS_B;
}   

double fast_opt_kernel_15_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return tanh(((*covw)*(*covw)*dot_prod_sparse_opt(x,y))+uc[0]);

    KILL_WARNINGS_B;
}

double fast_opt_kernel_16_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return tanh(((*covw)*(*covw)*dot_prod_sparse_opt(x,y))+uv[0]);

    KILL_WARNINGS_B;
}   

double fast_opt_kernel_17_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return exp((*covw)*(*covw)*dot_prod_sparse_opt(x,y));

    KILL_WARNINGS_B;
}

double fast_opt_kernel_18_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return exp((*covw)*(*covw)*dot_prod_sparse_opt(x,y))-uc[0];

    KILL_WARNINGS_B;
}

double fast_opt_kernel_19_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return exp((*covw)*(*covw)*dot_prod_sparse_opt(x,y))-uv[0];

    KILL_WARNINGS_B;
}

double fast_opt_kernel_20_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return exp(-((*covw)*(*covw)*dif_meas_sparse_opt(x,y)));

    KILL_WARNINGS_B;
}

double fast_opt_kernel_21_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return exp(-((*covw)*(*covw)*dif_meas_sparse_opt(x,y)))-uc[0];

    KILL_WARNINGS_B;
}

double fast_opt_kernel_22_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return exp(-((*covw)*(*covw)*dif_meas_sparse_opt(x,y)))-uv[0];

    KILL_WARNINGS_B;
}

double fast_opt_kernel_23_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return exp(-((*covw)*eucldist_sparse_opt(x,y)));

    KILL_WARNINGS_B;
}

double fast_opt_kernel_24_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return exp(-((*covw)*eucldist_sparse_opt(x,y)))-uc[0];

    KILL_WARNINGS_B;
}

double fast_opt_kernel_25_sparse(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw)
{
    return exp(-((*covw)*eucldist_sparse_opt(x,y)))-uv[0];

    KILL_WARNINGS_B;
}


opt_kern_ptr kernel_dicto::get_fast_kern(L_DOUBLE **uc, L_DOUBLE **uv, double **covw)
{
    opt_kern_ptr result = NULL;

    if ( dsize == 1 )
    {
        if ( kweights[0] == 1.0 )
        {
            if ( ( result = (kern[0])->optim_kernel_fast ) != NULL )
            {
                *uv   = &(((kern[0])->uv)[0]);
                *uc   = &(((kern[0])->uc)[0]);
                *covw = (&(kern[0])->covw);
            }
        }
    }

    return result;
}
