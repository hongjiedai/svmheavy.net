
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
// Generic SV data class
//
// Written by: Alistair Shilton
//             Melbourne University
//

#ifndef _svdata_h
#define _svdata_h

#include "svflags.h"
#include "kcache.h"
#include "svdata.h"
#include "kernel.h"
#include "common/vector.h"
#include "common/matrix.h"
#include "common/factor.h"

class SVdata;
class SVdata_short_summary;
class SVdata_summary;

std::ostream &operator<<(std::ostream &output, SVdata &dumpee);
std::istream &operator>>(std::istream &input,  SVdata &source);

class SVdata
{
    friend std::ostream &operator<<(std::ostream &, SVdata &);
    friend std::istream &operator>>(std::istream &, SVdata &);

    public:

    /*
       Constructors:

       SVdata(f,...)     - construct a variable bias dataset, linear kernel.
       SVdata(K,f,...)   - construct a variable bias dataset.
       SVdata(K,b,f,...) - construct a fixed bias dataset.

       in each case, f is the flags argument (see svflags.h) and ... is
       the triple cache_memsize, cache_min_rowdin, d2c_buff_size.

       cache_memsize    = the amount of memory to be allocated to the dynamic 
                          kernel cache, if it is used (in MB).
       cache_min_rowdim = minimum row size used when calculating the number
                          of rows stored in the dynamic kernel cache, if it
                          is used.
       d2c_buff_size    = size of buffer used by d2c solver.
    */

    SVdata(                                        int svflags, long _cache_memsize, long _cache_min_rowdim, long d2c_buff_size);
    SVdata(const kernel_dicto &_K,                 int svflags, long _cache_memsize, long _cache_min_rowdim, long d2c_buff_size);
    SVdata(const kernel_dicto &_K, L_DOUBLE b_fix, int svflags, long _cache_memsize, long _cache_min_rowdim, long d2c_buff_size);
    SVdata(const SVdata &source);
    SVdata(const SVdata_summary &_source);

    SVdata &operator=(const SVdata &right_op);

    L_DOUBLE calc_g(const fVECTOR &y);
    L_DOUBLE calc_g_grad(const fVECTOR &y);

    /*
       Manipulate training set.
    */

    void addpoint(int contype_i,
                  long ident_i,
                  const fVECTOR  &x_i,
                  const L_DOUBLE &z_i,
                  const L_DOUBLE &rho_i,
                  const L_DOUBLE &rho_star_i,
                  const L_DOUBLE &h_i,
                  const L_DOUBLE &v_i,
                  const L_DOUBLE &gamma_i,
                  const L_DOUBLE &gamma_star_i,
                  const L_DOUBLE &mu_i,
                  const L_DOUBLE &mu_star_i);
    void delpoint(long i);
    void force_zero(long i);

    void set_fixed_bias(const L_DOUBLE &new_bias);
    void set_variable_bias(void);
    void set_kernel(const kernel_dicto &_K);

    void set_z(long i, const L_DOUBLE &zval);
    void set_rho(long i, const L_DOUBLE &rhoval);
    void set_rho_star(long i, const L_DOUBLE &rho_starval);
    void set_h(long i, const L_DOUBLE &hval);
    void set_v(long i, const L_DOUBLE &vval);
    void set_gamma(long i, const L_DOUBLE &gammaval);
    void set_gamma_star(long i, const L_DOUBLE &gamma_starval);
    void set_mu(long i, const L_DOUBLE &muval);
    void set_mu_star(long i, const L_DOUBLE &mu_starval);

    void scale_z(const L_DOUBLE &zscale);
    void scale_rho(const L_DOUBLE &rhoscale);
    void scale_rho_star(const L_DOUBLE &rho_starscale);
    void scale_rho_both(const L_DOUBLE &rho_bothscale);
    void scale_h(const L_DOUBLE &hscale);
    void scale_v(const L_DOUBLE &vscale);
    void scale_hv(const L_DOUBLE &vhscale);
    void scale_gamma(const L_DOUBLE &gammascale);
    void scale_gamma_star(const L_DOUBLE &gamma_starscale);
    void scale_gamma_both(const L_DOUBLE &gamma_bothscale);
    void scale_mu(const L_DOUBLE &muscale);
    void scale_mu_star(const L_DOUBLE &mu_starscale);
    void scale_mu_both(const L_DOUBLE &mu_bothscale);

    /*
       Stuff.
    */

    L_DOUBLE get_alpha_pivot(long i_pivot);

    L_DOUBLE get_e(long i);

    void step_none(L_DOUBLE d_b);
    void step_one(L_DOUBLE d_alpha, long i);
    void step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore);
    void step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b);

    L_DOUBLE get_G_no_extras(long i, long j);
    L_DOUBLE get_G_tau(long i, long j);

    void constrain_Z(long i);
    void constrain_L(long i);
    void constrain_U(long i);
    void free_L(long i);
    void free_U(long i);
    void constrain_Z_pivot(long i_pivot);
    void constrain_L_pivot(long i_pivot);
    void constrain_U_pivot(long i_pivot);
    void free_L_pivot(long i_pivot);
    void free_U_pivot(long i_pivot);

    long find_marker(long id);
    int is_sparse(void) const;

    /*
       Optimisation functions.

       enter_opt must be called before using any of these functions, as
       all are reliant on fast pointers, which can be invalidated through
       the addition/deletion of points over time.

       exit_opt must be called after the optimisation is finished, to make
       any cumulative adjustments which have been postponed during
       optimisation to make things as fast as possible.

       at entry: initialise _fast vector references.
       at exit:  fix length of e_free_pivot vector (if unused) to match N_F.
                 if pivotting was off during optimisation, fix order vector.
    */

    void enter_opt(void);
    void exit_opt(void);

    L_DOUBLE opt_get_e(long iz);
    L_DOUBLE opt_get_e_pos(long iz);
    L_DOUBLE opt_get_e_neg(long iz);

    L_DOUBLE opt_get_e_pivot(long i_pivotz);
    L_DOUBLE opt_get_e_pos_pivot(long i_pivotz);
    L_DOUBLE opt_get_e_neg_pivot(long i_pivotz);

    fVECTOR &opt_get_e_free_pivot(void);

    L_DOUBLE opt_get_alpha_pivot(long i_pivotz);
    L_DOUBLE opt_get_h_pivot(long i_pivotz);
    L_DOUBLE opt_get_v_pivot(long i_pivotz);

    int opt_get_tau_pivot(long i_pivotz);
    int opt_get_contype_pivot(long i_pivotz);

    L_DOUBLE opt_get_G_tau(long iz, long jz);
    L_DOUBLE opt_get_G_extras(long iz, long jz);
    L_DOUBLE opt_get_G_tau_pivot(long i_pivotz, long j_pivotz);
    L_DOUBLE opt_calc_G_tau_from_scratch(long iz, long jz);

    void opt_step_none(L_DOUBLE d_b);
    void opt_step_one(L_DOUBLE d_alpha, long iz);
    void opt_step_two(L_DOUBLE d_b, L_DOUBLE d_f, L_DOUBLE d_alpha1, L_DOUBLE d_alpha2, long i1z, long i2z);
    void opt_step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore);
    void opt_step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b);

    void opt_fix_e_free_pivot(void);

    /*
       Pivotting.
    */

    void opt_constrain_Z(long iz);
    void opt_constrain_L(long iz);
    void opt_constrain_U(long iz);
    void opt_free_L(long iz);
    void opt_free_U(long iz);
    void opt_constrain_Z_pivot(long i_pivotz);
    void opt_constrain_L_pivot(long i_pivotz);
    void opt_constrain_U_pivot(long i_pivotz);
    void opt_free_L_pivot(long i_pivotz);
    void opt_free_U_pivot(long i_pivotz);

    /*
       Factorisation operations
    */

    long opt_get_nbad_factor(void);
    void opt_minverse_factor(fVECTOR &d_alpha_pivot, L_DOUBLE &d_b, long z_start, long z_end, int f_zero);
    void opt_minverse_factor(fVECTOR &d_alpha_pivot,                long z_start, long z_end            );
    void opt_near_invert_factor(fVECTOR &d_alpha_pivot, L_DOUBLE &d_b);
    void opt_near_invert_factor(fVECTOR &d_alpha_pivot               );

    /*
       Summaries
    */

    void regen_from_short_summary(const SVdata_short_summary &source);
    void remake_from_summary(const SVdata_summary &source);

    /*
       NB: data is nominally read-only from hereon.
    */

    /*
       Abstract training set.
    */

    f_fVECTOR x;
    fVECTOR z;
    fVECTOR rho;
    fVECTOR rho_star;
    fVECTOR v;
    fVECTOR h;
    fVECTOR gamma;
    fVECTOR gamma_star;
    fVECTOR mu;
    fVECTOR mu_star;
    int fix_bias;
    kernel_dicto K;

    /*
       Miscellaneous stuff.

       contype = 0: restrain to 0
                 1: restrain negative
                 2: restrain positive
                 3: unrestrained
    */

    iVECTOR contype;
    lVECTOR ident;
    fVECTOR free_ones;

    /*
       Abstract solution set.

       tau = 0:  constrained at 0
             +1: free positive
             -1: free negative
             +2: constrained at h
             -2: constrained at -v
    */

    fVECTOR alpha;
    L_DOUBLE b;
    iVECTOR tau;
    lVECTOR m;
    lVECTOR minv;

    /*
       Optimisation data.

       e_tau is the not quite the same as in the thesis.  In particular,
       if contype[i] == 1,2 and e_tau[i] is cached, the relevant rho
       factor will be added at all times, regardless of tau[i].
    */

    L_DOUBLE E;
    fVECTOR e_tau;
    L_DOUBLE f;
    fMATRIX G_tau;
    fMATRFact beth_tau;

    /*
       Fast pointers
    */

    L_DOUBLE  *fast_z;
    L_DOUBLE  *fast_rho;
    L_DOUBLE  *fast_rho_star;
    L_DOUBLE  *fast_v;
    L_DOUBLE  *fast_h;
    L_DOUBLE  *fast_gamma;
    L_DOUBLE  *fast_gamma_star;
    L_DOUBLE  *fast_mu;
    L_DOUBLE  *fast_mu_star;
    L_DOUBLE  *fast_alpha;
    L_DOUBLE  *fast_e_tau;
    int       *fast_contype;
    int       *fast_tau;
    long      *fast_m;
    long      *fast_minv;
    fVECTOR  **fast_x;
    L_DOUBLE **fast_opt_x;

    /*
       Ordered optimisation copy data
    */

    fVECTOR e_free_pivot;

    /*
       Problem size and stuff
    */

    long N,N_Z,N_L,N_U,N_F,N_FP,N_FN,N_C,N_S;
    int svflags;

    void fix_e_free_pivot(void);
    L_DOUBLE calc_G_tau_from_scratch(long i, long j);
    void rewrite_G_tau(int pad_first);

    /*
       Kernel caching stuff.
    */

    Kcache kerncache;
    opt_kern_ptr fast_kernel_fn;
    L_DOUBLE *fast_uv;
    L_DOUBLE *fast_uc;
    L_DOUBLE *fast_covw;

    long cache_memsize;
    long cache_min_rowdim;
    long d2c_buff_size;

    /*
       Function pointers.
    */

    L_DOUBLE (SVdata::*opt_sel__get_e)(long iz);
    L_DOUBLE (SVdata::*opt_sel__get_e_pos)(long iz);
    L_DOUBLE (SVdata::*opt_sel__get_e_neg)(long iz);
    L_DOUBLE (SVdata::*opt_sel__get_G_tau)(long iz, long jz);
    L_DOUBLE (SVdata::*opt_sel__get_G_extras)(long iz, long jz);
    L_DOUBLE (SVdata::*opt_sel__get_G_tau_pivot)(long i_pivotz, long j_pivotz);
    L_DOUBLE (SVdata::*opt_sel__calc_G_tau_from_scratch)(long iz, long jz);
    void     (SVdata::*opt_sel__step_none)(L_DOUBLE d_b);
    void     (SVdata::*opt_sel__step_one)(L_DOUBLE d_alpha, long iz);
    void     (SVdata::*opt_sel__step_two)(L_DOUBLE d_b, L_DOUBLE d_f, L_DOUBLE d_alpha1, L_DOUBLE d_alpha2, long i1z, long i2z);
    void     (SVdata::*opt_sel__step_general_pivot)(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore);
    void     (SVdata::*opt_sel__step_general_pivotx)(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b);
    void     (SVdata::*opt_sel__constrain_Z)(long iz);
    void     (SVdata::*opt_sel__constrain_L)(long iz);
    void     (SVdata::*opt_sel__constrain_U)(long iz);
    void     (SVdata::*opt_sel__free_L)(long iz);
    void     (SVdata::*opt_sel__free_U)(long iz);
    void     (SVdata::*opt_sel__constrain_Z_pivot)(long i_pivotz);
    void     (SVdata::*opt_sel__constrain_L_pivot)(long i_pivotz);
    void     (SVdata::*opt_sel__constrain_U_pivot)(long i_pivotz);
    void     (SVdata::*opt_sel__free_L_pivot)(long i_pivotz);
    void     (SVdata::*opt_sel__free_U_pivot)(long i_pivotz);

    /*
       Optimisation options.

       NB: indices of the form iz, i2z jnumz etc begin at 0 (not 1)

       Trivial functions: opt_d2csmo_get_G_extras
                          opt_d2csmo_constrain_L
                          opt_d2csmo_constrain_U
                          opt_d2csmo_free_L
                          opt_d2csmo_free_U
                          opt_d2csmo_pattern_get_e_pos
                          opt_d2csmo_pattern_get_e_neg
                          opt_d2csmo_pattern_get_G_extras
                          opt_d2csmo_pattern_constrain_Z
                          opt_d2csmo_pattern_constrain_L
                          opt_d2csmo_pattern_constrain_U
                          opt_d2csmo_pattern_free_L
                          opt_d2csmo_pattern_free_U
                          opt_d2csmo_pattern_fastkern_get_e_pos
                          opt_d2csmo_pattern_fastkern_get_e_neg
                          opt_d2csmo_pattern_fastkern_kfast_get_G_tau
                          opt_d2csmo_pattern_fastkern_get_G_extras
                          opt_d2csmo_pattern_fastkern_constrain_Z
                          opt_d2csmo_pattern_fastkern_constrain_L
                          opt_d2csmo_pattern_fastkern_constrain_U
                          opt_d2csmo_pattern_fastkern_free_L
                          opt_d2csmo_pattern_fastkern_free_U
                          opt_actsmall_get_e
                          opt_actsmall_get_G_tau
                          opt_actsmall_get_G_extras
                          opt_actmedium_get_e
                          opt_actmedium_get_G_extras

       Nonworking functions: opt_d2csmo_get_G_tau_pivot
                             opt_d2csmo_step_general_pivot
                             opt_d2csmo_step_general_pivotx
                             opt_d2csmo_constrain_Z_pivot
                             opt_d2csmo_constrain_L_pivot
                             opt_d2csmo_constrain_U_pivot
                             opt_d2csmo_free_L_pivot
                             opt_d2csmo_free_U_pivot
                             opt_d2csmo_pattern_get_G_tau_pivot
                             opt_d2csmo_pattern_step_general_pivot
                             opt_d2csmo_pattern_step_general_pivotx
                             opt_d2csmo_pattern_constrain_Z_pivot
                             opt_d2csmo_pattern_constrain_L_pivot
                             opt_d2csmo_pattern_constrain_U_pivot
                             opt_d2csmo_pattern_free_L_pivot
                             opt_d2csmo_pattern_free_U_pivot
                             opt_actsmall_step_none
                             opt_actsmall_step_one
                             opt_actsmall_step_two
                             opt_actsmall_constrain_Z
                             opt_actsmall_constrain_L
                             opt_actsmall_constrain_U
                             opt_actsmall_free_L
                             opt_actsmall_free_U
                             opt_actmedium_step_none
                             opt_actmedium_step_one
                             opt_actmedium_step_two
                             opt_actmedium_constrain_Z
                             opt_actmedium_constrain_L
                             opt_actmedium_constrain_U
                             opt_actmedium_free_L
                             opt_actmedium_free_U
       
    */

    L_DOUBLE opt_generic_get_e(long iz);
    L_DOUBLE opt_generic_get_e_pos(long iz);
    L_DOUBLE opt_generic_get_e_neg(long iz);
    L_DOUBLE opt_generic_get_G_tau(long iz, long jz);
    L_DOUBLE opt_generic_get_G_extras(long iz, long jz);
    L_DOUBLE opt_generic_get_G_tau_pivot(long i_pivotz, long j_pivotz);
    L_DOUBLE opt_generic_calc_G_tau_from_scratch(long iz, long jz);
    void     opt_generic_step_none(L_DOUBLE d_b);
    void     opt_generic_step_one(L_DOUBLE d_alpha, long iz);
    void     opt_generic_step_two(L_DOUBLE d_b, L_DOUBLE d_f, L_DOUBLE d_alpha1, L_DOUBLE d_alpha2, long i1z, long i2z);
    void     opt_generic_step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore);
    void     opt_generic_step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b);
    void     opt_generic_constrain_Z(long iz);
    void     opt_generic_constrain_L(long iz);
    void     opt_generic_constrain_U(long iz);
    void     opt_generic_free_L(long iz);
    void     opt_generic_free_U(long iz);
    void     opt_generic_constrain_Z_pivot(long i_pivotz);
    void     opt_generic_constrain_L_pivot(long i_pivotz);
    void     opt_generic_constrain_U_pivot(long i_pivotz);
    void     opt_generic_free_L_pivot(long i_pivotz);
    void     opt_generic_free_U_pivot(long i_pivotz);

    L_DOUBLE opt_actsmall_get_e(long iz);
    L_DOUBLE opt_actsmall_get_e_pos(long iz);
    L_DOUBLE opt_actsmall_get_e_neg(long iz);
    L_DOUBLE opt_actsmall_get_G_tau(long iz, long jz);
    L_DOUBLE opt_actsmall_get_G_extras(long iz, long jz);
    L_DOUBLE opt_actsmall_get_G_tau_pivot(long i_pivotz, long j_pivotz);
    L_DOUBLE opt_actsmall_calc_G_tau_from_scratch(long iz, long jz);
    void     opt_actsmall_step_none(L_DOUBLE d_b);
    void     opt_actsmall_step_one(L_DOUBLE d_alpha, long iz);
    void     opt_actsmall_step_two(L_DOUBLE d_b, L_DOUBLE d_f, L_DOUBLE d_alpha1, L_DOUBLE d_alpha2, long i1z, long i2z);
    void     opt_actsmall_step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore);
    void     opt_actsmall_step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b);
    void     opt_actsmall_constrain_Z(long iz);
    void     opt_actsmall_constrain_L(long iz);
    void     opt_actsmall_constrain_U(long iz);
    void     opt_actsmall_free_L(long iz);
    void     opt_actsmall_free_U(long iz);
    void     opt_actsmall_constrain_Z_pivot(long i_pivotz);
    void     opt_actsmall_constrain_L_pivot(long i_pivotz);
    void     opt_actsmall_constrain_U_pivot(long i_pivotz);
    void     opt_actsmall_free_L_pivot(long i_pivotz);
    void     opt_actsmall_free_U_pivot(long i_pivotz);

    L_DOUBLE opt_actmedium_get_e(long iz);
    L_DOUBLE opt_actmedium_get_e_pos(long iz);
    L_DOUBLE opt_actmedium_get_e_neg(long iz);
    L_DOUBLE opt_actmedium_get_G_tau(long iz, long jz);
    L_DOUBLE opt_actmedium_get_G_extras(long iz, long jz);
    L_DOUBLE opt_actmedium_get_G_tau_pivot(long i_pivotz, long j_pivotz);
    L_DOUBLE opt_actmedium_calc_G_tau_from_scratch(long iz, long jz);
    void     opt_actmedium_step_none(L_DOUBLE d_b);
    void     opt_actmedium_step_one(L_DOUBLE d_alpha, long iz);
    void     opt_actmedium_step_two(L_DOUBLE d_b, L_DOUBLE d_f, L_DOUBLE d_alpha1, L_DOUBLE d_alpha2, long i1z, long i2z);
    void     opt_actmedium_step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore);
    void     opt_actmedium_step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b);
    void     opt_actmedium_constrain_Z(long iz);
    void     opt_actmedium_constrain_L(long iz);
    void     opt_actmedium_constrain_U(long iz);
    void     opt_actmedium_free_L(long iz);
    void     opt_actmedium_free_U(long iz);
    void     opt_actmedium_constrain_Z_pivot(long i_pivotz);
    void     opt_actmedium_constrain_L_pivot(long i_pivotz);
    void     opt_actmedium_constrain_U_pivot(long i_pivotz);
    void     opt_actmedium_free_L_pivot(long i_pivotz);
    void     opt_actmedium_free_U_pivot(long i_pivotz);

    L_DOUBLE opt_d2csmo_get_e(long iz);
    L_DOUBLE opt_d2csmo_get_e_pos(long iz);
    L_DOUBLE opt_d2csmo_get_e_neg(long iz);
    L_DOUBLE opt_d2csmo_get_G_tau(long iz, long jz);
    L_DOUBLE opt_d2csmo_get_G_extras(long iz, long jz);
    L_DOUBLE opt_d2csmo_get_G_tau_pivot(long i_pivotz, long j_pivotz);
    L_DOUBLE opt_d2csmo_calc_G_tau_from_scratch(long iz, long jz);
    void     opt_d2csmo_step_none(L_DOUBLE d_b);
    void     opt_d2csmo_step_one(L_DOUBLE d_alpha, long iz);
    void     opt_d2csmo_step_two(L_DOUBLE d_b, L_DOUBLE d_f, L_DOUBLE d_alpha1, L_DOUBLE d_alpha2, long i1z, long i2z);
    void     opt_d2csmo_step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore);
    void     opt_d2csmo_step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b);
    void     opt_d2csmo_constrain_Z(long iz);
    void     opt_d2csmo_constrain_L(long iz);
    void     opt_d2csmo_constrain_U(long iz);
    void     opt_d2csmo_free_L(long iz);
    void     opt_d2csmo_free_U(long iz);
    void     opt_d2csmo_constrain_Z_pivot(long i_pivotz);
    void     opt_d2csmo_constrain_L_pivot(long i_pivotz);
    void     opt_d2csmo_constrain_U_pivot(long i_pivotz);
    void     opt_d2csmo_free_L_pivot(long i_pivotz);
    void     opt_d2csmo_free_U_pivot(long i_pivotz);

    L_DOUBLE opt_d2csmo_pattern_get_e(long iz);
    L_DOUBLE opt_d2csmo_pattern_get_e_pos(long iz);
    L_DOUBLE opt_d2csmo_pattern_get_e_neg(long iz);
    L_DOUBLE opt_d2csmo_pattern_get_G_tau(long iz, long jz);
    L_DOUBLE opt_d2csmo_pattern_get_G_extras(long iz, long jz);
    L_DOUBLE opt_d2csmo_pattern_get_G_tau_pivot(long i_pivotz, long j_pivotz);
    L_DOUBLE opt_d2csmo_pattern_calc_G_tau_from_scratch(long iz, long jz);
    void     opt_d2csmo_pattern_step_none(L_DOUBLE d_b);
    void     opt_d2csmo_pattern_step_one(L_DOUBLE d_alpha, long iz);
    void     opt_d2csmo_pattern_step_two(L_DOUBLE d_b, L_DOUBLE d_f, L_DOUBLE d_alpha1, L_DOUBLE d_alpha2, long i1z, long i2z);
    void     opt_d2csmo_pattern_step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore);
    void     opt_d2csmo_pattern_step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b);
    void     opt_d2csmo_pattern_constrain_Z(long iz);
    void     opt_d2csmo_pattern_constrain_L(long iz);
    void     opt_d2csmo_pattern_constrain_U(long iz);
    void     opt_d2csmo_pattern_free_L(long iz);
    void     opt_d2csmo_pattern_free_U(long iz);
    void     opt_d2csmo_pattern_constrain_Z_pivot(long i_pivotz);
    void     opt_d2csmo_pattern_constrain_L_pivot(long i_pivotz);
    void     opt_d2csmo_pattern_constrain_U_pivot(long i_pivotz);
    void     opt_d2csmo_pattern_free_L_pivot(long i_pivotz);
    void     opt_d2csmo_pattern_free_U_pivot(long i_pivotz);

    L_DOUBLE opt_d2csmo_pattern_kfast_get_e(long iz);
    L_DOUBLE opt_d2csmo_pattern_kfast_get_e_pos(long iz);
    L_DOUBLE opt_d2csmo_pattern_kfast_get_e_neg(long iz);
    L_DOUBLE opt_d2csmo_pattern_kfast_get_G_tau(long iz, long jz);
    L_DOUBLE opt_d2csmo_pattern_kfast_get_G_extras(long iz, long jz);
    L_DOUBLE opt_d2csmo_pattern_kfast_get_G_tau_pivot(long i_pivotz, long j_pivotz);
    L_DOUBLE opt_d2csmo_pattern_kfast_calc_G_tau_from_scratch(long iz, long jz);
    void     opt_d2csmo_pattern_kfast_step_none(L_DOUBLE d_b);
    void     opt_d2csmo_pattern_kfast_step_one(L_DOUBLE d_alpha, long iz);
    void     opt_d2csmo_pattern_kfast_step_two(L_DOUBLE d_b, L_DOUBLE d_f, L_DOUBLE d_alpha1, L_DOUBLE d_alpha2, long i1z, long i2z);
    void     opt_d2csmo_pattern_kfast_step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore);
    void     opt_d2csmo_pattern_kfast_step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b);
    void     opt_d2csmo_pattern_kfast_constrain_Z(long iz);
    void     opt_d2csmo_pattern_kfast_constrain_L(long iz);
    void     opt_d2csmo_pattern_kfast_constrain_U(long iz);
    void     opt_d2csmo_pattern_kfast_free_L(long iz);
    void     opt_d2csmo_pattern_kfast_free_U(long iz);
    void     opt_d2csmo_pattern_kfast_constrain_Z_pivot(long i_pivotz);
    void     opt_d2csmo_pattern_kfast_constrain_L_pivot(long i_pivotz);
    void     opt_d2csmo_pattern_kfast_constrain_U_pivot(long i_pivotz);
    void     opt_d2csmo_pattern_kfast_free_L_pivot(long i_pivotz);
    void     opt_d2csmo_pattern_kfast_free_U_pivot(long i_pivotz);
};

std::ostream &operator<<(std::ostream &output, SVdata_summary &dumpee);
std::istream &operator>>(std::istream &input,  SVdata_summary &source);

class SVdata_summary
{
    friend std::ostream &operator<<(std::ostream &, SVdata_summary &);
    friend std::istream &operator>>(std::istream &, SVdata_summary &);

    public:

    SVdata_summary(const SVdata &source);

    f_fVECTOR x;
    fVECTOR z;
    fVECTOR rho;
    fVECTOR rho_star;
    fVECTOR v;
    fVECTOR h;
    fVECTOR gamma;
    fVECTOR gamma_star;
    fVECTOR mu;
    fVECTOR mu_star;
    int fix_bias;
    kernel_dicto K;
    lVECTOR ident;
    fVECTOR free_ones;

    iVECTOR contype;
    fVECTOR alpha;
    L_DOUBLE b;
    iVECTOR tau;
    lVECTOR m;
    lVECTOR minv;
    L_DOUBLE E;
    fVECTOR e_tau;
    L_DOUBLE f;
    fMATRFact beth_tau;
    fVECTOR e_free_pivot;
    long N,N_Z,N_L,N_U,N_F,N_FP,N_FN,N_C,N_S;
    int svflags;

    long cache_memsize;
    long cache_min_rowdim;
    long d2c_buff_size;
};


std::ostream &operator<<(std::ostream &output, SVdata_short_summary &dumpee);

class SVdata_short_summary
{
    friend std::ostream &operator<<(std::ostream &, SVdata_short_summary &);

    public:

    SVdata_short_summary(const SVdata &source);

    iVECTOR contype;
    fVECTOR alpha;
    L_DOUBLE b;
    iVECTOR tau;
    lVECTOR m;
    lVECTOR minv;
    L_DOUBLE E;
    fVECTOR e_tau;
    L_DOUBLE f;
    fMATRFact beth_tau;
    fVECTOR e_free_pivot;
    long N,N_Z,N_L,N_U,N_F,N_FP,N_FN,N_C,N_S;
    int svflags;

    long cache_memsize;
    long cache_min_rowdim;
    long d2c_buff_size;
};

#endif

