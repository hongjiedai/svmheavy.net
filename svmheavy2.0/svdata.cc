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


#include "common/outfilt.h"

/*
   All operations including beth_tau (the hessian factorisation) which are
   themselves explicitly dependent on G_tau have been placed inside the
   following macros.
*/

#define FIX_FACT_BETH                                                   \
{                                                                       \
    if ( fix_bias )                                                     \
    {                                                                   \
        if ( KERN_CACHE_FULL(svflags) )                                 \
        {                                                               \
            if ( H_FACT_CHOL(svflags) )                                 \
            {                                                           \
                m.set_offsets(N_C+1,N);                                 \
                fMATRFact temp(1,m,G_tau,ZZERO,FACT_CHOL);              \
                m.reset_offsets();                                      \
                                                                        \
                beth_tau = temp;                                        \
            }                                                           \
                                                                        \
            else if ( H_FACT_INVE(svflags) )                            \
            {                                                           \
                m.set_offsets(N_C+1,N);                                 \
                fMATRFact temp(1,m,G_tau,ZZERO,FACT_INVERSE);           \
                m.reset_offsets();                                      \
                                                                        \
                beth_tau = temp;                                        \
            }                                                           \
        }                                                               \
                                                                        \
        else if ( KERN_CACHE_FREE(svflags) )                            \
        {                                                               \
            if ( H_FACT_CHOL(svflags) )                                 \
            {                                                           \
                fMATRFact temp(0,m,G_tau,ZZERO,FACT_CHOL);              \
                                                                        \
                beth_tau = temp;                                        \
            }                                                           \
                                                                        \
            else if ( H_FACT_INVE(svflags) )                            \
            {                                                           \
                fMATRFact temp(0,m,G_tau,ZZERO,FACT_INVERSE);           \
                                                                        \
                beth_tau = temp;                                        \
            }                                                           \
        }                                                               \
    }                                                                   \
                                                                        \
    else                                                                \
    {                                                                   \
        if ( KERN_CACHE_FULL(svflags) )                                 \
        {                                                               \
            if ( H_FACT_CHOL(svflags) )                                 \
            {                                                           \
                m.set_offsets(N_C+1,N);                                 \
                free_ones.set_offsets(N_C+1,N);                         \
                fMATRFact temp(1,m,G_tau,free_ones,ZZERO,FACT_CHOL);    \
                free_ones.reset_offsets();                              \
                m.reset_offsets();                                      \
                                                                        \
                beth_tau = temp;                                        \
            }                                                           \
                                                                        \
            else if ( H_FACT_INVE(svflags) )                            \
            {                                                           \
                m.set_offsets(N_C+1,N);                                 \
                free_ones.set_offsets(N_C+1,N);                         \
                fMATRFact temp(1,m,G_tau,free_ones,ZZERO,FACT_INVERSE); \
                free_ones.reset_offsets();                              \
                m.reset_offsets();                                      \
                                                                        \
                beth_tau = temp;                                        \
            }                                                           \
        }                                                               \
                                                                        \
        else if ( KERN_CACHE_FREE(svflags) )                            \
        {                                                               \
            if ( H_FACT_CHOL(svflags) )                                 \
            {                                                           \
                free_ones.set_offsets(N_C+1,N);                         \
                fMATRFact temp(0,m,G_tau,free_ones,ZZERO,FACT_CHOL);    \
                free_ones.reset_offsets();                              \
                                                                        \
                beth_tau = temp;                                        \
            }                                                           \
                                                                        \
            else if ( H_FACT_INVE(svflags) )                            \
            {                                                           \
                free_ones.set_offsets(N_C+1,N);                         \
                fMATRFact temp(0,m,G_tau,free_ones,ZZERO,FACT_INVERSE); \
                free_ones.reset_offsets();                              \
                                                                        \
                beth_tau = temp;                                        \
            }                                                           \
        }                                                               \
    }                                                                   \
}

#define PERT_ONE_BETH(__i,__pert)                                       \
{                                                                       \
    if ( fix_bias )                                                     \
    {                                                                   \
        m.set_offsets(N_C+1,N);                                         \
        beth_tau.pert_one(m,__i,__pert,G_tau);                          \
        m.reset_offsets();                                              \
    }                                                                   \
                                                                        \
    else                                                                \
    {                                                                   \
        m.set_offsets(N_C+1,N);                                         \
        free_ones.set_offsets(N_C+1,N);                                 \
        beth_tau.pert_one(m,__i,__pert,G_tau,free_ones);                \
        free_ones.reset_offsets();                                      \
        m.reset_offsets();                                              \
    }                                                                   \
}

#define RANKONE_BETH(__vect,__scale)                                    \
{                                                                       \
    if ( fix_bias )                                                     \
    {                                                                   \
        m.set_offsets(N_C+1,N);                                         \
        beth_tau.rankone(m,__vect,__scale,G_tau);                       \
        m.reset_offsets();                                              \
    }                                                                   \
                                                                        \
    else                                                                \
    {                                                                   \
        m.set_offsets(N_C+1,N);                                         \
        free_ones.set_offsets(N_C+1,N);                                 \
        beth_tau.rankone(m,__vect,__scale,G_tau,free_ones);             \
        free_ones.reset_offsets();                                      \
        m.reset_offsets();                                              \
    }                                                                   \
}

#define SHRINK_BETH(__i)                                                \
{                                                                       \
    if ( fix_bias )                                                     \
    {                                                                   \
        m.set_offsets(N_C+1,N);                                         \
        beth_tau.fact_shrink(m,__i,G_tau);                              \
        m.reset_offsets();                                              \
    }                                                                   \
                                                                        \
    else                                                                \
    {                                                                   \
        m.set_offsets(N_C+1,N);                                         \
        free_ones.set_offsets(N_C+1,N);                                 \
        beth_tau.fact_shrink(m,__i,G_tau,free_ones);                    \
        free_ones.reset_offsets();                                      \
        m.reset_offsets();                                              \
    }                                                                   \
}

#define ADDEND_BETH(__vect)                                             \
{                                                                       \
    if ( fix_bias )                                                     \
    {                                                                   \
        m.set_offsets(N_C+1,N);                                         \
        beth_tau.fact_addend(m,__vect,G_tau);                           \
        m.reset_offsets();                                              \
    }                                                                   \
                                                                        \
    else                                                                \
    {                                                                   \
        m.set_offsets(N_C+1,N);                                         \
        free_ones.set_offsets(N_C+1,N);                                 \
        beth_tau.fact_addend(m,__vect,1.0,G_tau,free_ones);             \
        free_ones.reset_offsets();                                      \
        m.reset_offsets();                                              \
    }                                                                   \
}

#define MINVERSE_BETH(__d_a_piv,__d_b,__e_piv,__f,__z_start,__z_end,__f_zero) \
{                                                                       \
    m.set_offsets(N_C+1,N);                                             \
    free_ones.set_offsets(N_C+1,N);                                     \
    beth_tau.minverse(m,__d_a_piv,__d_b,__e_piv,__f,G_tau,free_ones,__z_start,__z_end,__f_zero); \
    free_ones.reset_offsets();                                          \
    m.reset_offsets();                                                  \
}

#define MINVERSE_BETH_SHORT(__d_a_piv,__e_piv,__z_start,__z_end)        \
{                                                                       \
    m.set_offsets(N_C+1,N);                                             \
    beth_tau.minverse(m,__d_a_piv,__e_piv,G_tau,__z_start,__z_end);     \
    m.reset_offsets();                                                  \
}

#define NEARINVERSE_BETH(__d_a_piv,__d_b)                               \
{                                                                       \
    m.set_offsets(N_C+1,N);                                             \
    free_ones.set_offsets(N_C+1,N);                                     \
    beth_tau.near_invert(m,__d_a_piv,__d_b,G_tau,free_ones);            \
    free_ones.reset_offsets();                                          \
    m.reset_offsets();                                                  \
}

#define NEARINVERSE_BETH_SHORT(__d_a_piv)                               \
{                                                                       \
    m.set_offsets(N_C+1,N);                                             \
    beth_tau.near_invert(m,__d_a_piv,G_tau);                            \
    m.reset_offsets();                                                  \
}





//
// Generic SV data class
//
// Written by: Alistair Shilton
//             Melbourne University
//

#include "svdata.h"

SVdata::SVdata(int _svflags, long _cache_memsize, long _cache_min_rowdim, long _d2c_buff_size)
{
    FN_ENTRY_POINT;

    fast_x     = NULL;
    fast_opt_x = NULL;

    N   = 0;
    N_Z = 0;
    N_L = 0;
    N_U = 0;
    N_F = 0;
    N_FP = 0;
    N_FN = 0;
    N_C = 0;
    N_S = 0;

    svflags = _svflags;

    contype.make_zero(DO_FORCE);
    z.make_zero(DO_FORCE);
    rho.make_zero(DO_FORCE);
    rho_star.make_zero(DO_FORCE);
    v.make_zero(DO_FORCE);
    h.make_zero(DO_FORCE);
    gamma.make_zero(DO_FORCE);
    gamma_star.make_zero(DO_FORCE);
    mu.make_zero(DO_FORCE);
    mu_star.make_zero(DO_FORCE);
    free_ones.make_zero(DO_FORCE);
    ident.make_zero(DO_FORCE),
    alpha.make_zero(DO_FORCE);
    tau.make_zero(DO_FORCE);
    m.make_zero(DO_FORCE);
    minv.make_zero(DO_FORCE);
    e_tau.make_zero(DO_FORCE);
    e_free_pivot.make_zero(DO_FORCE);

    contype.make_normal(DO_FORCE);
    z.make_normal(DO_FORCE);
    rho.make_normal(DO_FORCE);
    rho_star.make_normal(DO_FORCE);
    v.make_normal(DO_FORCE);
    h.make_normal(DO_FORCE);
    gamma.make_normal(DO_FORCE);
    gamma_star.make_normal(DO_FORCE);
    mu.make_normal(DO_FORCE);
    mu_star.make_normal(DO_FORCE);
    free_ones.make_normal(DO_FORCE);
    ident.make_normal(DO_FORCE),
    alpha.make_normal(DO_FORCE);
    tau.make_normal(DO_FORCE);
    m.make_normal(DO_FORCE);
    minv.make_normal(DO_FORCE);
    e_tau.make_normal(DO_FORCE);
    e_free_pivot.make_normal(DO_FORCE);

    fix_bias = 0;

    b = 0.0;
    E = 0.0;
    f = 0.0;

    G_tau.make_symmetric(DO_FORCE);

    FIX_FACT_BETH;

    cache_memsize    = _cache_memsize;
    cache_min_rowdim = _cache_min_rowdim;

    d2c_buff_size = _d2c_buff_size;

    if ( KERN_CACHE_DYNA(svflags) )
    {
        kerncache.setmemsize(this,cache_memsize,cache_min_rowdim);
    }

    rewrite_G_tau(1);

    FN_EXIT_POINT;
}

SVdata::SVdata(const kernel_dicto &_K, int _svflags, long _cache_memsize, long _cache_min_rowdim, long _d2c_buff_size)
{
    FN_ENTRY_POINT;

    fast_x     = NULL;
    fast_opt_x = NULL;

    N   = 0;
    N_Z = 0;
    N_L = 0;
    N_U = 0;
    N_F = 0;
    N_FP = 0;
    N_FN = 0;
    N_C = 0;
    N_S = 0;

    svflags = _svflags;

    K = _K;

    contype.make_zero(DO_FORCE);
    z.make_zero(DO_FORCE);
    rho.make_zero(DO_FORCE);
    rho_star.make_zero(DO_FORCE);
    v.make_zero(DO_FORCE);
    h.make_zero(DO_FORCE);
    gamma.make_zero(DO_FORCE);
    gamma_star.make_zero(DO_FORCE);
    mu.make_zero(DO_FORCE);
    mu_star.make_zero(DO_FORCE);
    free_ones.make_zero(DO_FORCE);
    ident.make_zero(DO_FORCE),
    alpha.make_zero(DO_FORCE);
    tau.make_zero(DO_FORCE);
    m.make_zero(DO_FORCE);
    minv.make_zero(DO_FORCE);
    e_tau.make_zero(DO_FORCE);
    e_free_pivot.make_zero(DO_FORCE);

    contype.make_normal(DO_FORCE);
    z.make_normal(DO_FORCE);
    rho.make_normal(DO_FORCE);
    rho_star.make_normal(DO_FORCE);
    v.make_normal(DO_FORCE);
    h.make_normal(DO_FORCE);
    gamma.make_normal(DO_FORCE);
    gamma_star.make_normal(DO_FORCE);
    mu.make_normal(DO_FORCE);
    mu_star.make_normal(DO_FORCE);
    free_ones.make_normal(DO_FORCE);
    ident.make_normal(DO_FORCE),
    alpha.make_normal(DO_FORCE);
    tau.make_normal(DO_FORCE);
    m.make_normal(DO_FORCE);
    minv.make_normal(DO_FORCE);
    e_tau.make_normal(DO_FORCE);
    e_free_pivot.make_normal(DO_FORCE);

    fix_bias = 0;

    b = 0.0;
    E = 0.0;
    f = 0.0;

    G_tau.make_symmetric(DO_FORCE);

    FIX_FACT_BETH;

    cache_memsize    = _cache_memsize;
    cache_min_rowdim = _cache_min_rowdim;

    d2c_buff_size = _d2c_buff_size;

    if ( KERN_CACHE_DYNA(svflags) )
    {
        kerncache.setmemsize(this,cache_memsize,cache_min_rowdim);
    }

    rewrite_G_tau(1);

    FN_EXIT_POINT;
}

SVdata::SVdata(const kernel_dicto &_K, L_DOUBLE b_fix, int _svflags, long _cache_memsize, long _cache_min_rowdim, long _d2c_buff_size)
{
    FN_ENTRY_POINT;

    fast_x     = NULL;
    fast_opt_x = NULL;

    N   = 0;
    N_Z = 0;
    N_L = 0;
    N_U = 0;
    N_F = 0;
    N_FP = 0;
    N_FN = 0;
    N_C = 0;
    N_S = 0;

    svflags = _svflags;

    K = _K;

    contype.make_zero(DO_FORCE);
    z.make_zero(DO_FORCE);
    rho.make_zero(DO_FORCE);
    rho_star.make_zero(DO_FORCE);
    v.make_zero(DO_FORCE);
    h.make_zero(DO_FORCE);
    gamma.make_zero(DO_FORCE);
    gamma_star.make_zero(DO_FORCE);
    mu.make_zero(DO_FORCE);
    mu_star.make_zero(DO_FORCE);
    free_ones.make_zero(DO_FORCE);
    ident.make_zero(DO_FORCE),
    alpha.make_zero(DO_FORCE);
    tau.make_zero(DO_FORCE);
    m.make_zero(DO_FORCE);
    minv.make_zero(DO_FORCE);
    e_tau.make_zero(DO_FORCE);
    e_free_pivot.make_zero(DO_FORCE);

    contype.make_normal(DO_FORCE);
    z.make_normal(DO_FORCE);
    rho.make_normal(DO_FORCE);
    rho_star.make_normal(DO_FORCE);
    v.make_normal(DO_FORCE);
    h.make_normal(DO_FORCE);
    gamma.make_normal(DO_FORCE);
    gamma_star.make_normal(DO_FORCE);
    mu.make_normal(DO_FORCE);
    mu_star.make_normal(DO_FORCE);
    free_ones.make_normal(DO_FORCE);
    ident.make_normal(DO_FORCE),
    alpha.make_normal(DO_FORCE);
    tau.make_normal(DO_FORCE);
    m.make_normal(DO_FORCE);
    minv.make_normal(DO_FORCE);
    e_tau.make_normal(DO_FORCE);
    e_free_pivot.make_normal(DO_FORCE);

    fix_bias = 1;

    b = b_fix;
    E = 0.0;
    f = 0.0;

    G_tau.make_symmetric(DO_FORCE);

    FIX_FACT_BETH;

    cache_memsize    = _cache_memsize;
    cache_min_rowdim = _cache_min_rowdim;

    d2c_buff_size = _d2c_buff_size;

    if ( KERN_CACHE_DYNA(svflags) )
    {
        kerncache.setmemsize(this,cache_memsize,cache_min_rowdim);
    }

    rewrite_G_tau(1);

    FN_EXIT_POINT;
}

SVdata::SVdata(const SVdata &source)
{
    FN_ENTRY_POINT;

    fast_x     = NULL;
    fast_opt_x = NULL;

    long i;

    fix_bias = source.fix_bias;
    svflags  = source.svflags;

    N    = source.N;
    N_Z  = source.N_Z;
    N_L  = source.N_L;
    N_U  = source.N_U;
    N_F  = source.N_F;
    N_FP = source.N_FP;
    N_FN = source.N_FN;
    N_C  = source.N_C;
    N_S  = source.N_S;

    b = source.b;
    E = source.E;
    f = source.f;

    x.trim_to_size(0);

    if ( N > 0 )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            x.addend((source.x).get_offset_element(i-1));
        }
    }

    z.make_zero(DO_FORCE);             z            = source.z;
    rho.make_zero(DO_FORCE);           rho          = source.rho;
    rho_star.make_zero(DO_FORCE);      rho_star     = source.rho_star;
    v.make_zero(DO_FORCE);             v            = source.v;
    h.make_zero(DO_FORCE);             h            = source.h;
    gamma.make_zero(DO_FORCE);         gamma        = source.gamma;
    gamma_star.make_zero(DO_FORCE);    gamma_star   = source.gamma_star;
    mu.make_zero(DO_FORCE);            mu           = source.mu;
    mu_star.make_zero(DO_FORCE);       mu_star      = source.mu_star;
    contype.make_zero(DO_FORCE);       contype      = source.contype;
    ident.make_zero(DO_FORCE);         ident        = source.ident;
    free_ones.make_zero(DO_FORCE);     free_ones    = source.free_ones;
    alpha.make_zero(DO_FORCE);         alpha        = source.alpha;
    tau.make_zero(DO_FORCE);           tau          = source.tau;
    m.make_zero(DO_FORCE);             m            = source.m;
    minv.make_zero(DO_FORCE);          minv         = source.minv;
    e_tau.make_zero(DO_FORCE);         e_tau        = source.e_tau;
    G_tau.make_zero(DO_FORCE);         G_tau        = source.G_tau;
                                       beth_tau     = source.beth_tau;
    e_free_pivot.make_zero(DO_FORCE);  e_free_pivot = source.e_free_pivot;

    K = source.K;

    cache_memsize    = source.cache_memsize;
    cache_min_rowdim = source.cache_min_rowdim;

    d2c_buff_size = source.d2c_buff_size;

    if ( KERN_CACHE_DYNA(svflags) )
    {
        kerncache.setmemsize(this,cache_memsize,cache_min_rowdim);
    }

    FN_EXIT_POINT;
}

SVdata::SVdata(const SVdata_summary &source)
{
    FN_ENTRY_POINT;

    fast_x     = NULL;
    fast_opt_x = NULL;

    long i;

    fix_bias = source.fix_bias;
    svflags  = source.svflags;

    N    = source.N;
    N_Z  = source.N_Z;
    N_L  = source.N_L;
    N_U  = source.N_U;
    N_F  = source.N_F;
    N_FP = source.N_FP;
    N_FN = source.N_FN;
    N_C  = source.N_C;
    N_S  = source.N_S;

    b = source.b;
    E = source.E;
    f = source.f;

    x.trim_to_size(0);

    if ( N > 0 )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            x.addend((source.x).get_offset_element(i-1));
        }
    }

    z.make_zero(DO_FORCE);             z            = source.z;
    rho.make_zero(DO_FORCE);           rho          = source.rho;
    rho_star.make_zero(DO_FORCE);      rho_star     = source.rho_star;
    v.make_zero(DO_FORCE);             v            = source.v;
    h.make_zero(DO_FORCE);             h            = source.h;
    gamma.make_zero(DO_FORCE);         gamma        = source.gamma;
    gamma_star.make_zero(DO_FORCE);    gamma_star   = source.gamma_star;
    mu.make_zero(DO_FORCE);            mu           = source.mu;
    mu_star.make_zero(DO_FORCE);       mu_star      = source.mu_star;
    contype.make_zero(DO_FORCE);       contype      = source.contype;
    ident.make_zero(DO_FORCE);         ident        = source.ident;
    free_ones.make_zero(DO_FORCE);     free_ones    = source.free_ones;
    alpha.make_zero(DO_FORCE);         alpha        = source.alpha;
    tau.make_zero(DO_FORCE);           tau          = source.tau;
    m.make_zero(DO_FORCE);             m            = source.m;
    minv.make_zero(DO_FORCE);          minv         = source.minv;
    e_tau.make_zero(DO_FORCE);         e_tau        = source.e_tau;
                                       beth_tau     = source.beth_tau;
    e_free_pivot.make_zero(DO_FORCE);  e_free_pivot = source.e_free_pivot;

    K = source.K;

    cache_memsize    = source.cache_memsize;
    cache_min_rowdim = source.cache_min_rowdim;

    d2c_buff_size = source.d2c_buff_size;

    if ( KERN_CACHE_DYNA(svflags) )
    {
        kerncache.setmemsize(this,cache_memsize,cache_min_rowdim);
    }

    /*
       Remake G_tau
    */

    rewrite_G_tau(1);

    FN_EXIT_POINT;
}

SVdata &SVdata::operator=(const SVdata &source)
{
    FN_ENTRY_POINT;

    long i;

    fix_bias = source.fix_bias;
    svflags  = source.svflags;

    N    = source.N;
    N_Z  = source.N_Z;
    N_L  = source.N_L;
    N_U  = source.N_U;
    N_F  = source.N_F;
    N_FP = source.N_FP;
    N_FN = source.N_FN;
    N_C  = source.N_C;
    N_S  = source.N_S;

    b = source.b;
    E = source.E;
    f = source.f;

    x.trim_to_size(0);

    if ( N > 0 )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            x.addend((source.x).get_offset_element(i-1));
        }
    }

                                       x            = source.x;
    z.make_zero(DO_FORCE);             z            = source.z;
    rho.make_zero(DO_FORCE);           rho          = source.rho;
    rho_star.make_zero(DO_FORCE);      rho_star     = source.rho_star;
    v.make_zero(DO_FORCE);             v            = source.v;
    h.make_zero(DO_FORCE);             h            = source.h;
    gamma.make_zero(DO_FORCE);         gamma        = source.gamma;
    gamma_star.make_zero(DO_FORCE);    gamma_star   = source.gamma_star;
    mu.make_zero(DO_FORCE);            mu           = source.mu;
    mu_star.make_zero(DO_FORCE);       mu_star      = source.mu_star;
    contype.make_zero(DO_FORCE);       contype      = source.contype;
    ident.make_zero(DO_FORCE);         ident        = source.ident;
    free_ones.make_zero(DO_FORCE);     free_ones    = source.free_ones;
    alpha.make_zero(DO_FORCE);         alpha        = source.alpha;
    tau.make_zero(DO_FORCE);           tau          = source.tau;
    m.make_zero(DO_FORCE);             m            = source.m;
    minv.make_zero(DO_FORCE);          minv         = source.minv;
    e_tau.make_zero(DO_FORCE);         e_tau        = source.e_tau;
    G_tau.make_zero(DO_FORCE);         G_tau        = source.G_tau;
                                       beth_tau     = source.beth_tau;
    e_free_pivot.make_zero(DO_FORCE);  e_free_pivot = source.e_free_pivot;

    K = source.K;

    cache_memsize    = source.cache_memsize;
    cache_min_rowdim = source.cache_min_rowdim;

    d2c_buff_size = source.d2c_buff_size;

    if ( KERN_CACHE_DYNA(svflags) )
    {
        kerncache.setmemsize(this,cache_memsize,cache_min_rowdim);
    }

    FN_EXIT_POINT (*this);
}

L_DOUBLE SVdata::calc_g(const fVECTOR &y)
{
    FN_ENTRY_POINT;

    L_DOUBLE result;
    long j_pivot;

    result = b;

    if ( N_Z < N )
    {
        for ( j_pivot = N_Z ; j_pivot < N ; j_pivot++ )
        {
            result += alpha[m[j_pivot]-1] * K.kernel(x[m[j_pivot]-1],y);
        }
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::calc_g_grad(const fVECTOR &y)
{
    FN_ENTRY_POINT;

    L_DOUBLE result;
    long j_pivot;

    result = 0.0;

    if ( N_Z < N )
    {
        for ( j_pivot = N_Z ; j_pivot < N ; j_pivot++ )
        {
            result += alpha[m[j_pivot]-1] * K.cross_kernel(x[m[j_pivot]-1],y);
        }
    }

    FN_EXIT_POINT result;
}

void SVdata::addpoint(int contype_i,
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
                     const L_DOUBLE &mu_star_i)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( contype_i >= 0 );
    THROW_ASSERT( contype_i <= 3 );
    THROW_ASSERT( v_i <= 0.0 );
    THROW_ASSERT( h_i >= 0.0 );
    THROW_ASSERT( gamma_i      >= 0.0 );
    THROW_ASSERT( gamma_star_i >= 0.0 );
    THROW_ASSERT( mu_i      >= 0.0 );
    THROW_ASSERT( mu_star_i >= 0.0 );
    THROW_ASSERT( ( contype_i != 3 ) || ( rho_i      >= 0.0 ) );
    THROW_ASSERT( ( contype_i != 3 ) || ( rho_star_i >= 0.0 ) );

    if ( ( rho_i != 0.0 ) || ( rho_star_i != 0.0 ) )
    {
        MAKE_SVM_HINT_RHO_NONZERO(svflags);
    }

    if ( ( gamma_i != 0.0 ) || ( gamma_star_i != 0.0 ) )
    {
        MAKE_SVM_HINT_GAMMA_NONZERO(svflags);
    }

    if ( ( mu_i != 0.0 ) || ( mu_star_i != 0.0 ) )
    {
        MAKE_SVM_HINT_MU_NONZERO(svflags);
    }

    long j,j_pivot;
    L_DOUBLE new_e;

    /*
       Add point to training data
    */

    x.addend(x_i);
    z.addend(z_i);
    rho.addend(rho_i);
    rho_star.addend(rho_star_i);
    v.addend(v_i);
    h.addend(h_i);
    gamma.addend(gamma_i);
    gamma_star.addend(gamma_star_i);
    mu.addend(mu_i);
    mu_star.addend(mu_star_i);

    contype.addend(contype_i);
    ident.addend(ident_i);
    free_ones.addend(1.0);

    /*
       When the point is added, alpha is always set to 0 (and therefore
       tau is zero also).
    */

    alpha.addend(0.0);
    tau.addend(0);

    /*
       Extend G_tau (if a complete kernel cache is kept)
    */

    if ( KERN_CACHE_FULL(svflags) )
    {
        fVECTOR Gend('x',N+1);

        for ( j = 0 ; j <= N ; j++ )
        {
            Gend[j] = K.kernel(x[N],x[j],NULL,NULL,0.0,0.0,N,j+1,z_i,z[j]);
        }

        G_tau.addend(Gend[N],Gend);
    }

    else if ( KERN_CACHE_DYNA(svflags) )
    {
        N++;
        N_Z++;
        N_C++;

        kerncache.setmemsize(this,cache_memsize,cache_min_rowdim);

        N--;
        N_Z--;
        N_C--;
    }

    /*
       Calculate e and extend e_tau
    */

    e_tau.addend(0.0);

    if ( E_CACHE_FULL(svflags) )
    {
        new_e = b-z_i;

        if ( N_Z < N )
        {
            for ( j_pivot = N_Z ; j_pivot < N ; j_pivot++ )
            {
                new_e += get_G_tau(N+1,m[j_pivot])*alpha[m[j_pivot]-1];
            }
        }

        if ( contype_i == 2 )
        {
            new_e += rho_i;
        }

        else if ( contype_i == 1 )
        {
            new_e -= rho_star_i;
        }

        e_tau[N] = new_e;
    }

    /*
       This vector is added to the START of the pivotted vectors, and hence
       we must add the pointer to the start of the order vector.
    */

    m.addstart(N+1);

    N++;
    N_Z++;
    N_C++;

    minv.addend(0);

    for ( j = 0 ; j < N ; j++ )
    {
        minv[j]++;
    }

    FN_EXIT_POINT;
}

void SVdata::delpoint(long i)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );

    long j_pivot;

    force_zero(i);

    /*
       Point has already been forced to zero, so need to remove it from
       G_tau only if a full cache of kernels is kept.
    */

    if ( KERN_CACHE_FULL(svflags) )
    {
        G_tau.remove(i);
    }

    else if ( KERN_CACHE_DYNA(svflags) )
    {
        kerncache.remove(this,i);
    }

    x.remove(i);
    z.remove(i);
    rho.remove(i);
    rho_star.remove(i);
    v.remove(i);
    h.remove(i);
    gamma.remove(i);
    gamma_star.remove(i);
    mu.remove(i);
    mu_star.remove(i);

    contype.remove(i);
    ident.remove(i);
    free_ones.remove(i);

    alpha.remove(i);
    tau.remove(i);
    e_tau.remove(i);

    /*
       Finally, the order vector will need to be re-numbered at the same
       time as the relevant element is removed.
    */

    for ( j_pivot = 0 ; j_pivot < N ; j_pivot++ )
    {
        if ( m[j_pivot] == i )
        {
            m.remove(j_pivot+1);

            j_pivot--;
            N--;
            N_Z--;
            N_C--;
        }

        else if ( m[j_pivot] > i )
        {
            m[j_pivot]--;
        }
    }

    minv.remove(1);

    if ( N > 0 )
    {
        for ( j_pivot = 0 ; j_pivot < N ; j_pivot++ )
        {
            minv[m[j_pivot]-1] = j_pivot+1;
        }
    }

    FN_EXIT_POINT;
}

void SVdata::force_zero(long i)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );

    /*
       Two stages: 1. if its at an upper/lower bound, free it.
                   2. then, if it's free, change alpha and constrain it to 0
    */

    if ( tau[i-1] == -2 )
    {
        /*
           Currently the point is constrained at a lower bound.  Hence we
           need to free it to a free negative.
        */

        free_L(i);
    }

    else if ( tau[i-1] == +2 )
    {
        /*
           Currently the point is constrained at an upper bound.  Hence we
           need to free it to a free positive.
        */

        free_U(i);
    }

    if ( tau[i-1] != 0 )
    {
        step_one(-alpha[i-1],i);
        constrain_Z(i);
    }

    /*
       Finally, set contype
    */

    contype[i-1] = 0;

    FN_EXIT_POINT;
}

void SVdata::set_fixed_bias(const L_DOUBLE &new_bias)
{
    FN_ENTRY_POINT;

    step_none(new_bias-b);

    fix_bias = 1;

    FIX_FACT_BETH;

    FN_EXIT_POINT;
}

void SVdata::set_variable_bias(void)
{
    FN_ENTRY_POINT;

    fix_bias = 0;

    FIX_FACT_BETH;

    FN_EXIT_POINT;
}

void SVdata::set_kernel(const kernel_dicto &_K)
{
    FN_ENTRY_POINT;

    long i_pivot,j_pivot;

    K = _K;

    /*
       fix G_tau - important to remember matrix symmetry here (don't want
       to add anything twice by mistake).
    */

    rewrite_G_tau(0);

    if ( N > 0 )
    {
        /*
           fix e_tau
        */

        L_DOUBLE new_e;
        long i_pivot_start;

        if ( E_CACHE_FULL(svflags) )
        {
            i_pivot_start = 1;
        }

        else if ( E_CACHE_SUPP(svflags) )
        {
            i_pivot_start = N_Z+1;
        }

        else if ( E_CACHE_FREE(svflags) )
        {
            i_pivot_start = N_C+1;
        }

        else
        {
            i_pivot_start = N+1;
        }

        if ( i_pivot_start <= N )
        {
            for ( i_pivot = i_pivot_start-1 ; i_pivot < N ; i_pivot++ )
            {
                new_e = b - z[m[i_pivot]-1];

                if ( N_Z < N )
                {
                    for ( j_pivot = N_Z ; j_pivot < N ; j_pivot++ )
                    {
                        new_e += get_G_tau(m[i_pivot],m[j_pivot]) * alpha[m[j_pivot]-1];
                    }
                }

                if ( ( tau[m[i_pivot]-1] > 0 ) || ( contype[m[i_pivot]-1] == 2 ) )
                {
                    new_e += rho[m[i_pivot]-1];
                }

                else if ( ( tau[m[i_pivot]-1] < 0 ) || ( contype[m[i_pivot]-1] == 1 ) )
                {
                    new_e -= rho_star[m[i_pivot]-1];
                }

                e_tau[m[i_pivot]-1] = new_e;
            }

            fix_e_free_pivot();
        }

        /*
           fix beth_tau
        */

        FIX_FACT_BETH;
    }

    FN_EXIT_POINT;
}

void SVdata::set_z(long i, const L_DOUBLE &zval)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );

    if ( E_CACHE_FULL(svflags) )
    {
        e_tau[i-1] += z[i-1];
        z[i-1] = zval;
        e_tau[i-1] -= z[i-1];

        if ( ( tau[i-1] == +1 ) || ( tau[i-1] == -1 ) )
        {
            fix_e_free_pivot();
        }
    }

    else if ( E_CACHE_SUPP(svflags) )
    {
        if ( tau[i-1] != 0 )
        {
            e_tau[i-1] += z[i-1];
            z[i-1] = zval;
            e_tau[i-1] -= z[i-1];

            if ( ( tau[i-1] == +1 ) || ( tau[i-1] == -1 ) )
            {
                fix_e_free_pivot();
            }
        }

        else
        {
            z[i-1] = zval;
        }
    }

    else if ( E_CACHE_FREE(svflags) )
    {
        if ( ( tau[i-1] == +1 ) || ( tau[i-1] == -1 ) )
        {
            e_tau[i-1] += z[i-1];
            z[i-1] = zval;
            e_tau[i-1] -= z[i-1];

            fix_e_free_pivot();
        }

        else
        {
            z[i-1] = zval;
        }
    }

    else
    {
        z[i-1] = zval;
    }

    FN_EXIT_POINT;
}

void SVdata::set_rho(long i, const L_DOUBLE &rhoval)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );
    THROW_ASSERT( ( rhoval >= 0.0 ) || ( contype[i-1] != 3 ) );

    if ( ( rhoval != 0.0 ) || !SVM_HINT_RHO_ZERO(svflags) )
    {
        MAKE_SVM_HINT_RHO_NONZERO(svflags);

        if ( ( tau[i-1] > 0 ) || ( contype[i-1] == 2 ) )
        {
            if ( E_CACHE_FULL(svflags) )
            {
                e_tau[i-1] -= rho[i-1];
                rho[i-1] = rhoval;
                e_tau[i-1] += rho[i-1];

                if ( tau[i-1] == +1 )
                {
                    fix_e_free_pivot();
                }
            }

            else if ( E_CACHE_SUPP(svflags) )
            {
                if ( tau[i-1] != 0 )
                {
                    e_tau[i-1] -= rho[i-1];
                    rho[i-1] = rhoval;
                    e_tau[i-1] += rho[i-1];

                    if ( tau[i-1] == +1 )
                    {
                        fix_e_free_pivot();
                    }
                }

                else
                {
                    rho[i-1] = rhoval;
                }
            }

            else if ( E_CACHE_FREE(svflags) )
            {
                if ( tau[i-1] == +1 )
                {
                    e_tau[i-1] -= rho[i-1];
                    rho[i-1] = rhoval;
                    e_tau[i-1] += rho[i-1];

                    fix_e_free_pivot();
                }

                else
                {
                    rho[i-1] = rhoval;
                }
            }

            else
            {
                rho[i-1] = rhoval;
            }
        }

        else
        {
            rho[i-1] = rhoval;
        }
    }

    FN_EXIT_POINT;
}

void SVdata::set_rho_star(long i, const L_DOUBLE &rho_starval)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );
    THROW_ASSERT( ( rho_starval >= 0.0 ) || ( contype[i-1] != 3 ) );

    if ( ( rho_starval != 0.0 ) || !SVM_HINT_RHO_ZERO(svflags) )
    {
        MAKE_SVM_HINT_RHO_NONZERO(svflags);

        if ( ( tau[i-1] < 0 ) || ( contype[i-1] == 1 ) )
        {
            if ( E_CACHE_FULL(svflags) )
            {
                e_tau[i-1] += rho_star[i-1];
                rho_star[i-1] = rho_starval;
                e_tau[i-1] -= rho_star[i-1];

                if ( tau[i-1] == -1 )
                {
                    fix_e_free_pivot();
                }
            }

            else if ( E_CACHE_SUPP(svflags) )
            {
                if ( tau[i-1] != 0 )
                {
                    e_tau[i-1] += rho_star[i-1];
                    rho_star[i-1] = rho_starval;
                    e_tau[i-1] -= rho_star[i-1];

                    if ( tau[i-1] == -1 )
                    {
                        fix_e_free_pivot();
                    }
                }

                else
                {
                    rho_star[i-1] = rho_starval;
                }
            }

            else if ( E_CACHE_FREE(svflags) )
            {
                if ( tau[i-1] == -1 )
                {
                    e_tau[i-1] += rho_star[i-1];
                    rho_star[i-1] = rho_starval;
                    e_tau[i-1] -= rho_star[i-1];

                    fix_e_free_pivot();
                }

                else
                {
                    rho_star[i-1] = rho_starval;
                }
            }

            else
            {
                rho_star[i-1] = rho_starval;
            }
        }

        else
        {
            rho_star[i-1] = rho_starval;
        }
    }

    FN_EXIT_POINT;
}

void SVdata::set_h(long i, const L_DOUBLE &hval)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );
    THROW_ASSERT( hval >= 0.0 );

    if ( tau[i-1] == +1 )
    {
        /*
           Not at a bound, but some care is required.  In particular, we
           need to make sure that alpha doesn't cross its bound.
        */

        if ( alpha[i-1] <= hval )
        {
            /*
               We're safe, can just set the bound with no affect.
            */

            h[i-1] = hval;
        }

        else
        {
            /*
               Will hit bound - first move bound to alpha, then constrain,
               then squeeze the bound.
            */

            set_h(i,alpha[i-1]);
            constrain_U(i);
            set_h(i,hval);
        }
    }

    else if ( tau[i-1] == +2 )
    {
        step_one(hval-alpha[i-1],i);
        alpha[i-1] = hval;
        h[i-1]     = hval;
    }

    else
    {
        h[i-1] = hval;
    }

    FN_EXIT_POINT;
}

void SVdata::set_v(long i, const L_DOUBLE &vval)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );
    THROW_ASSERT( vval <= 0.0 );

    if ( tau[i-1] == -1 )
    {
        /*
           Not at a bound, but some care is required.  In particular, we
           need to make sure that alpha doesn't cross its bound.
        */

        if ( alpha[i-1] >= vval )
        {
            /*
               We're safe, can just set the bound with no affect.
            */

            v[i-1] = vval;
        }

        else
        {
            /*
               Will hit bound - first move bound to alpha, then constrain,
               then squeeze the bound.
            */

            set_v(i,alpha[i-1]);
            constrain_L(i);
            set_v(i,vval);
        }
    }

    else if ( tau[i-1] == -2 )
    {
        /*
           Adjusting the bound now will affect alpha, so need to adjust lots
           of other stuff as well.
        */

        step_one(vval-alpha[i-1],i);
        alpha[i-1] = vval;
        v[i-1]     = vval;
    }

    else
    {
        v[i-1] = vval;
    }

    FN_EXIT_POINT;
}

void SVdata::set_gamma(long i, const L_DOUBLE &gammaval)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );
    THROW_ASSERT( gammaval >= 0.0 );

    if ( ( gammaval != 0.0 ) || !SVM_HINT_GAMMA_ZERO(svflags) )
    {
        MAKE_SVM_HINT_GAMMA_NONZERO(svflags);

        if ( tau[i-1] > 0 )
        {
            if ( ( tau[i-1] == +1 ) && ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) ) )
            {
                /*
                   fix beth_tau
                */

                PERT_ONE_BETH(minv[i-1],-gamma[i-1]);
            }

            if ( ( E_CACHE_FREE(svflags) && ( tau[i-1] == +1 ) ) || E_CACHE_SUPP(svflags) || E_CACHE_FULL(svflags) )
            {
                e_tau[i-1] -= gamma[i-1] * alpha[i-1];
            }

            if ( KERN_CACHE_FULL(svflags) )
            {
                G_tau(i-1,i-1) -= gamma[i-1];
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                if ( tau[i+1] == +1 )
                {
                    G_tau(minv[i-1]-N_C-1,minv[i-1]-N_C-1) -= gamma[i-1];
                }
            }

            gamma[i-1] = gammaval;

            if ( ( tau[i-1] == +1 ) && ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) ) )
            {
                /*
                   fix beth_tau
                */

                PERT_ONE_BETH(minv[i-1],gamma[i-1]);
            }

            if ( KERN_CACHE_FULL(svflags) )
            {
                G_tau(i-1,i-1) += gamma[i-1];
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                if ( tau[i+1] == +1 )
                {
                    G_tau(minv[i-1]-N_C-1,minv[i-1]-N_C-1) += gamma[i-1];
                }
            }

            if ( ( E_CACHE_FREE(svflags) && ( tau[i-1] == +1 ) ) || E_CACHE_SUPP(svflags) || E_CACHE_FULL(svflags) )
            {
                e_tau[i-1] += gamma[i-1] * alpha[i-1];
    
                if ( tau[i-1] == +1 )
                {
                    fix_e_free_pivot();
                }
            }
        }

        else
        {
            gamma[i-1] = gammaval;
        }
    }

    FN_EXIT_POINT;
}

void SVdata::set_gamma_star(long i, const L_DOUBLE &gamma_starval)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );
    THROW_ASSERT( gamma_starval >= 0.0 );

    if ( ( gamma_starval != 0.0 ) || !SVM_HINT_GAMMA_ZERO(svflags) )
    {
        MAKE_SVM_HINT_GAMMA_NONZERO(svflags);

        if ( tau[i-1] < 0 )
        {
            if ( ( tau[i-1] == -1 ) && ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) ) )
            {
                /*
                   fix beth_tau
                */

                PERT_ONE_BETH(minv[i-1],-gamma_star[i-1]);
            }

            if ( ( E_CACHE_FREE(svflags) && ( tau[i-1] == -1 ) ) || E_CACHE_SUPP(svflags) || E_CACHE_FULL(svflags) )
            {
                e_tau[i-1] -= gamma_star[i-1] * alpha[i-1];
            }

            if ( KERN_CACHE_FULL(svflags) )
            {
                G_tau(i-1,i-1) -= gamma_star[i-1];
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                if ( tau[i-1] == -1 )
                {
                    G_tau(minv[i-1]-N_C-1,minv[i-1]-N_C-1) -= gamma_star[i-1];
                }
            }

            gamma_star[i-1] = gamma_starval;

            if ( ( tau[i-1] == -1 ) && ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) ) )
            {
                /*
                   fix beth_tau
                */

                PERT_ONE_BETH(minv[i-1],gamma_star[i-1]);
            }

            if ( KERN_CACHE_FULL(svflags) )
            {
                G_tau(i-1,i-1) += gamma_star[i-1];
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                if ( tau[i-1] == -1 )
                {
                    G_tau(minv[i-1]-N_C-1,minv[i-1]-N_C-1) += gamma_star[i-1];
                }
            }

            if ( ( E_CACHE_FREE(svflags) && ( tau[i-1] == -1 ) ) || E_CACHE_SUPP(svflags) || E_CACHE_FULL(svflags) )
            {
                e_tau[i-1] += gamma_star[i-1] * alpha[i-1];

                if ( tau[i-1] == -1 )
                {
                    fix_e_free_pivot();
                }
            }
        }

        else
        {
            gamma_star[i-1] = gamma_starval;
        }
    }

    FN_EXIT_POINT;
}

void SVdata::set_mu(long i, const L_DOUBLE &muval)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );
    THROW_ASSERT( muval >= 0.0 );

    long j;
    long j_pivot;
    long j_pivot_start;

    if ( ( muval != 0.0 ) || !SVM_HINT_MU_ZERO(svflags) )
    {
        MAKE_SVM_HINT_MU_NONZERO(svflags);

        if ( tau[i-1] > 0 )
        {
            /*
               There are 3 affected variables here, namely E, G_tau and e_tau.
               First, we subtract the mu affects from each (as if mu were
               being set to zero), then we add them back in with the correct
               values.
            */

            if ( ( ( tau[i-1] == +1 ) || ( tau[i-1] == -1 ) ) && ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) ) )
            {
                /*
                   fix beth_tau
                */

                fVECTOR mu_x('c',N_F);

                for ( j_pivot = N_C ; j_pivot < N ; j_pivot++ )
                {
                    if ( tau[m[j_pivot]-1] == +1 )
                    {
                        mu_x[j_pivot] = mu[m[j_pivot]-1];
                    }

                    else
                    {
                        mu_x[j_pivot] = -mu_star[m[j_pivot]-1];
                    }
                }

                RANKONE_BETH(mu_x,-1.0);
            }

            E -= mu[i-1] * alpha[i-1];

            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j = 0 ; j < N ; j++ )
                {
                    if ( tau[j] > 0 )
                    {
                        G_tau(i-1,j) -= mu[j] * mu[i-1];
                    }

                    else if ( tau[j] < 0 )
                    {
                        G_tau(i-1,j) += mu_star[j] * mu[i-1];
                    }
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                if ( tau[i+1] == +1 )
                {
                    for ( j_pivot = N_C ; j_pivot < N ; j_pivot++ )
                    {
                        if ( tau[m[j_pivot]-1] > 0 )
                        {
                            G_tau(minv[i-1]-N_C-1,j_pivot-N_C) -= mu[m[j_pivot]-1] * mu[i-1];
                        }

                        else if ( tau[m[j_pivot]-1] < 0 )
                        {
                            G_tau(minv[i-1]-N_C-1,j_pivot-N_C) += mu_star[m[j_pivot]-1] * mu[i-1];
                        }
                    }
                }
            }

            if ( E_CACHE_FULL(svflags) )
            {
                j_pivot_start = 1;
            }

            else if ( E_CACHE_SUPP(svflags) )
            {
                j_pivot_start = N_Z+1;
            }

            else if ( E_CACHE_FREE(svflags) )
            {
                j_pivot_start = N_C+1;
            }

            else
            {
                j_pivot_start = N+1;
            }

            if ( j_pivot_start <= N )
            {
                for ( j_pivot = j_pivot_start-1 ; j_pivot < N ; j_pivot++ )
                {
                    if ( tau[m[j_pivot]-1] > 0 )
                    {
                        e_tau[m[j_pivot]-1] -= mu[m[j_pivot]-1] * mu[i-1] * alpha[i-1];
                    }

                    else if ( tau[m[j_pivot-1]-1] < 0 )
                    {
                        e_tau[m[j_pivot]-1] += mu_star[m[j_pivot]-1] * mu[i-1] * alpha[i-1];
                    }
                }
            }

            mu[i-1] = muval;

            if ( ( ( tau[i-1] == +1 ) || ( tau[i-1] == -1 ) ) && ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) ) )
            {
                /*
                   fix beth_tau
                */

                fVECTOR mu_x('c',N_F);

                for ( j_pivot = N_C ; j_pivot < N ; j_pivot++ )
                {
                    if ( tau[m[j_pivot]-1] == +1 )
                    {
                        mu_x[j_pivot] = mu[m[j_pivot]-1];
                    }

                    else
                    {
                        mu_x[j_pivot] = -mu_star[m[j_pivot]-1];
                    }
                }

                RANKONE_BETH(mu_x,+1.0);
            }

            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j = 0 ; j < N ; j++ )
                {
                    if ( tau[j] > 0 )
                    {
                        G_tau(i-1,j) += mu[j] * mu[i-1];
                    }

                    else if ( tau[j-1] < 0 )
                    {
                        G_tau(i-1,j) -= mu_star[j] * mu[i-1];
                    }
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                if ( tau[i+1] == +1 )
                {
                    for ( j_pivot = N_C+1 ; j_pivot <= N ; j_pivot++ )
                    {
                        if ( tau[m[j_pivot-1]-1] > 0 )
                        {
                            G_tau(minv[i-1]-N_C-1,j_pivot-N_C-1) += mu[m[j_pivot-1]-1] * mu[i-1];
                        }

                        else if ( tau[m[j_pivot-1]-1] < 0 )
                        {
                            G_tau(minv[i-1]-N_C-1,j_pivot-N_C-1) -= mu_star[m[j_pivot-1]-1] * mu[i-1];
                        }
                    }
                }
            }

            if ( j_pivot_start <= N )
            {
                for ( j_pivot = j_pivot_start ; j_pivot <= N ; j_pivot++ )
                {
                    if ( tau[m[j_pivot-1]-1] > 0 )
                    {
                        e_tau[m[j_pivot-1]-1] += mu[m[j_pivot-1]-1] * mu[i-1] * alpha[i-1];
                    }

                    else if ( tau[m[j_pivot-1]-1] < 0 )
                    {
                        e_tau[m[j_pivot-1]-1] -= mu_star[m[j_pivot-1]-1] * mu[i-1] * alpha[i-1];
                    }
                }

                fix_e_free_pivot();
            }

            E += mu[i-1] * alpha[i-1];
        }

        else
        {
            mu[i-1] = muval;
        }
    }

    FN_EXIT_POINT;
}

void SVdata::set_mu_star(long i, const L_DOUBLE &mu_starval)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );
    THROW_ASSERT( mu_starval >= 0.0 );

    long j;
    long j_pivot;
    long j_pivot_start;

    if ( ( mu_starval != 0.0 ) || !SVM_HINT_MU_ZERO(svflags) )
    {
        MAKE_SVM_HINT_MU_NONZERO(svflags);

        if ( tau[i-1] < 0 )
        {
            /*
               There are 3 affected variables here, namely E, G_tau and e_tau.
               First, we subtract the mu affects from each (as if mu were
               being set to zero), then we add them back in with the correct
               values.
            */

            if ( ( ( tau[i-1] == +1 ) || ( tau[i-1] == -1 ) ) && ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) ) )
            {
                /*
                   fix beth_tau
                */

                fVECTOR mu_x('c',N_F);

                for ( j_pivot = N_C+1 ; j_pivot <= N ; j_pivot++ )
                {
                    if ( tau[m[j_pivot-1]-1] == +1 )
                    {
                        mu_x[j_pivot-1] = -mu[m[j_pivot-1]-1];
                    }

                    else
                    {
                        mu_x[j_pivot-1] = mu_star[m[j_pivot-1]-1];
                    }
                }

                RANKONE_BETH(mu_x,-1.0);
            }

            E += mu_star[i-1] * alpha[i-1];

            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j = 1 ; j <= N ; j++ )
                {
                    if ( tau[j-1] > 0 )
                    {
                        G_tau(i-1,j-1) += mu[j-1] * mu_star[i-1];
                    }

                    else if ( tau[j-1] < 0 )
                    {
                        G_tau(i-1,j-1) -= mu_star[j-1] * mu_star[i-1];
                    }
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                if ( tau[i+1] == -1 )
                {
                    for ( j_pivot = N_C+1 ; j_pivot <= N ; j_pivot++ )
                    {
                        if ( tau[m[j_pivot-1]-1] > 0 )
                        {
                            G_tau(minv[i-1]-N_C-1,j_pivot-N_C-1) += mu[m[j_pivot-1]-1] * mu_star[i-1];
                        }

                        else if ( tau[m[j_pivot-1]-1] < 0 )
                        {
                            G_tau(minv[i-1]-N_C-1,j_pivot-N_C-1) -= mu_star[m[j_pivot-1]-1] * mu_star[i-1];
                        }
                    }
                }
            }

            if ( E_CACHE_FULL(svflags) )
            {
                j_pivot_start = 1;
            }

            else if ( E_CACHE_SUPP(svflags) )
            {
                j_pivot_start = N_Z+1;
            }

            else if ( E_CACHE_FREE(svflags) )
            {
                j_pivot_start = N_C+1;
            }

            else
            {
                j_pivot_start = N+1;
            }

            if ( j_pivot_start <= N )
            {
                for ( j_pivot = j_pivot_start ; j_pivot <= N ; j_pivot++ )
                {
                    if ( tau[m[j_pivot-1]-1] > 0 )
                    {
                        e_tau[m[j_pivot-1]-1] += mu[m[j_pivot-1]-1] * mu_star[i-1] * alpha[i-1];
                    }

                    else if ( tau[m[j_pivot-1]-1] < 0 )
                    {
                        e_tau[m[j_pivot-1]-1] -= mu_star[m[j_pivot-1]-1] * mu_star[i-1] * alpha[i-1];
                    }
                }
            }

            mu_star[i-1] = mu_starval;

            if ( ( ( tau[i-1] == +1 ) || ( tau[i-1] == -1 ) ) && ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) ) )
            {
                /*
                   fix beth_tau
                */

                fVECTOR mu_x('c',N_F);

                for ( j_pivot = N_C+1 ; j_pivot <= N ; j_pivot++ )
                {
                    if ( tau[m[j_pivot-1]-1] == +1 )
                    {
                        mu_x[j_pivot-1] = +mu[m[j_pivot-1]-1];
                    }

                    else
                    {
                        mu_x[j_pivot-1] = mu_star[m[j_pivot-1]-1];
                    }
                }

                RANKONE_BETH(mu_x,+1.0);
            }

            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j = 1 ; j <= N ; j++ )
                {
                    if ( tau[j-1] > 0 )
                    {
                        G_tau(i-1,j-1) -= mu[j-1] * mu_star[i-1];
                    }

                    else if ( tau[j-1] < 0 )
                    {
                        G_tau(i-1,j-1) += mu_star[j-1] * mu_star[i-1];
                    }
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                if ( tau[i+1] == -1 )
                {
                    for ( j_pivot = N_C+1 ; j_pivot <= N ; j_pivot++ )
                    {
                        if ( tau[m[j_pivot-1]-1] > 0 )
                        {
                            G_tau(minv[i-1]-N_C-1,j_pivot-N_C-1) -= mu[m[j_pivot-1]-1] * mu_star[i-1];
                        }

                        else if ( tau[m[j_pivot-1]-1] < 0 )
                        {
                            G_tau(minv[i-1]-N_C-1,j_pivot-N_C-1) += mu_star[m[j_pivot-1]-1] * mu_star[i-1];
                        }
                    }
                }
            }

            if ( j_pivot_start <= N )
            {
                for ( j_pivot = j_pivot_start ; j_pivot <= N ; j_pivot++ )
                {
                    if ( tau[m[j_pivot-1]-1] > 0 )
                    {
                        e_tau[m[j_pivot-1]-1] -= mu[m[j_pivot-1]-1] * mu_star[i-1] * alpha[i-1];
                    }

                    else if ( tau[m[j_pivot-1]-1] < 0 )
                    {
                        e_tau[m[j_pivot-1]-1] += mu_star[m[j_pivot-1]-1] * mu_star[i-1] * alpha[i-1];
                    }
                }

                fix_e_free_pivot();
            }

            E -= mu_star[i-1] * alpha[i-1];
        }

        else
        {
            mu_star[i-1] = mu_starval;
        }
    }

    FN_EXIT_POINT;
}

void SVdata::scale_z(const L_DOUBLE &zscale)
{
    FN_ENTRY_POINT;

    if ( E_CACHE_FULL(svflags) )
    {
        e_tau -= z;
        z *= zscale;
        e_tau += z;

        fix_e_free_pivot();
    }

    else if ( E_CACHE_SUPP(svflags) )
    {
        /*
           Technically, only a subset of e_tau should be updated here.
           However, it's just too much effort for too little result.
        */

        e_tau -= z;
        z *= zscale;
        e_tau += z;

        fix_e_free_pivot();
    }

    else if ( E_CACHE_FREE(svflags) )
    {
        /*
           Technically, only a subset of e_tau should be updated here.
           However, it's just too much effort for too little result.
        */

        e_tau -= z;
        z *= zscale;
        e_tau += z;

        fix_e_free_pivot();
    }

    else
    {
        z *= zscale;
    }

    FN_EXIT_POINT;
}

void SVdata::scale_rho(const L_DOUBLE &rhoscale)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( rhoscale >= 0.0 );

    long i;

    if ( ( rhoscale != 1.0 ) && !SVM_HINT_RHO_ZERO(svflags) )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            set_rho(i,rhoscale*rho[i-1]);
        }
    }

    if ( rhoscale == 0.0 )
    {
        MAKE_SVM_HINT_RHO_ZERO(svflags);
    }

    FN_EXIT_POINT;
}

void SVdata::scale_rho_star(const L_DOUBLE &rho_starscale)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( rho_starscale >= 0.0 );

    long i;

    if ( ( rho_starscale != 1.0 ) && !SVM_HINT_RHO_ZERO(svflags) )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            set_rho_star(i,rho_starscale*rho_star[i-1]);
        }
    }

    if ( rho_starscale == 0.0 )
    {
        MAKE_SVM_HINT_RHO_ZERO(svflags);
    }

    FN_EXIT_POINT;
}

void SVdata::scale_rho_both(const L_DOUBLE &rho_bothscale)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( rho_bothscale >= 0.0 );

    long i;

    if ( ( rho_bothscale != 1.0 ) && !SVM_HINT_RHO_ZERO(svflags) )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            set_rho(i,rho_bothscale*rho[i-1]);
            set_rho_star(i,rho_bothscale*rho_star[i-1]);
        }
    }

    if ( rho_bothscale == 0.0 )
    {
        MAKE_SVM_HINT_RHO_ZERO(svflags);
    }

    FN_EXIT_POINT;
}

void SVdata::scale_h(const L_DOUBLE &hscale)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( hscale >= 0.0 );

    long i;

    if ( hscale != 1.0 )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            set_h(i,hscale*h[i-1]);
        }
    }

    FN_EXIT_POINT;
}

void SVdata::scale_v(const L_DOUBLE &vscale)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( vscale >= 0.0 );

    long i;

    if ( vscale != 1.0 )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            set_v(i,vscale*v[i-1]);
        }
    }

    FN_EXIT_POINT;
}

void SVdata::scale_hv(const L_DOUBLE &vhscale)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( vhscale >= 0.0 );

    long i;

    if ( vhscale != 1.0 )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            set_h(i,vhscale*h[i-1]);
            set_v(i,vhscale*v[i-1]);
        }
    }

    FN_EXIT_POINT;
}

void SVdata::scale_gamma(const L_DOUBLE &gammascale)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( gammascale >= 0.0 );

    long i;

    if ( ( gammascale != 1.0 ) && !SVM_HINT_GAMMA_ZERO(svflags) )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            set_gamma(i,gammascale*gamma[i-1]);
        }
    }

    if ( gammascale == 0.0 )
    {
        MAKE_SVM_HINT_GAMMA_ZERO(svflags);
    }

    FN_EXIT_POINT;
}

void SVdata::scale_gamma_star(const L_DOUBLE &gamma_starscale)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( gamma_starscale >= 0.0 );

    long i;

    if ( ( gamma_starscale != 1.0 ) && !SVM_HINT_GAMMA_ZERO(svflags) )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            set_gamma_star(i,gamma_starscale*gamma_star[i-1]);
        }
    }

    if ( gamma_starscale == 0.0 )
    {
        MAKE_SVM_HINT_GAMMA_ZERO(svflags);
    }

    FN_EXIT_POINT;
}

void SVdata::scale_gamma_both(const L_DOUBLE &gamma_bothscale)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( gamma_bothscale >= 0.0 );

    long i;

    if ( ( gamma_bothscale != 1.0 ) && !SVM_HINT_GAMMA_ZERO(svflags) )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            set_gamma(i,gamma_bothscale*gamma[i-1]);
            set_gamma_star(i,gamma_bothscale*gamma_star[i-1]);
        }
    }

    if ( gamma_bothscale == 0.0 )
    {
        MAKE_SVM_HINT_GAMMA_ZERO(svflags);
    }

    FN_EXIT_POINT;
}

void SVdata::scale_mu(const L_DOUBLE &muscale)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( muscale >= 0.0 );

    long i;

    if ( ( muscale != 1.0 ) && !SVM_HINT_MU_ZERO(svflags) )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            set_mu(i,muscale*mu[i-1]);
        }
    }

    if ( muscale == 0.0 )
    {
        MAKE_SVM_HINT_MU_ZERO(svflags);
    }

    FN_EXIT_POINT;
}

void SVdata::scale_mu_star(const L_DOUBLE &mu_starscale)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( mu_starscale >= 0.0 );

    long i;

    if ( ( mu_starscale != 1.0 ) && !SVM_HINT_MU_ZERO(svflags) )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            set_mu_star(i,mu_starscale*mu_star[i-1]);
        }
    }

    if ( mu_starscale == 0.0 )
    {
        MAKE_SVM_HINT_MU_ZERO(svflags);
    }

    FN_EXIT_POINT;
}

void SVdata::scale_mu_both(const L_DOUBLE &mu_bothscale)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( mu_bothscale >= 0.0 );

    long i;

    if ( ( mu_bothscale != 1.0 ) && !SVM_HINT_MU_ZERO(svflags) )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            set_mu(i,mu_bothscale*mu[i-1]);
            set_mu_star(i,mu_bothscale*mu_star[i-1]);
        }
    }

    if ( mu_bothscale == 0.0 )
    {
        MAKE_SVM_HINT_MU_ZERO(svflags);
    }

    FN_EXIT_POINT;
}



/*
   stuff
*/

L_DOUBLE SVdata::get_alpha_pivot(long i_pivot)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivot >= 1 );
    THROW_ASSERT( i_pivot <= N );

    FN_EXIT_POINT alpha[m[i_pivot-1]-1];
}

L_DOUBLE SVdata::get_e(long i)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );

    L_DOUBLE result;

    if ( E_CACHE_FULL(svflags) )
    {
        result = e_tau[i-1];

        if ( ( tau[i-1] == 0 ) && ( contype[i-1] == 2 ) )
        {
            e_tau[i-1] -= rho[i-1];
        }

        else if ( ( tau[i-1] == 0 ) && ( contype[i-1] == 1 ) )
        {
            e_tau[i-1] += rho_star[i-1];
        }
    }

    else if ( E_CACHE_SUPP(svflags) )
    {
        if ( tau[i-1] != 0 )
        {
            result = e_tau[i-1];
        }

        else
        {
            goto calc_e_from_scratch;
        }
    }

    else if ( E_CACHE_FREE(svflags) )
    {
        if ( ( tau[i-1] == +1 ) || ( tau[i-1] == -1 ) )
        {
            result = e_tau[i-1];
        }

        else
        {
            goto calc_e_from_scratch;
        }
    }

    else
    {
        calc_e_from_scratch:

        long j_pivot;

        result = b - z[i-1];

        if ( N_Z < N )
        {
            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j_pivot = N_Z ; j_pivot < N ; j_pivot++ )
                {
                    result += G_tau(i-1,m[j_pivot]-1) * alpha[m[j_pivot]-1];
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                long i_pivot;

                i_pivot = minv[i-1];

                if ( i_pivot <= N_C )
                {
                    for ( j_pivot = N_Z ; j_pivot < N ; j_pivot++ )
                    {
                        result += get_G_tau(i,m[j_pivot]) * alpha[m[j_pivot]-1];
                    }
                }

                else
                {
                    if ( N_Z < N_C )
                    {
                        for ( j_pivot = N_Z ; j_pivot < N_C ; j_pivot++ )
                        {
                            result += get_G_tau(i,m[j_pivot]) * alpha[m[j_pivot]-1];
                        }
                    }

                    for ( j_pivot = N_C ; j_pivot < N ; j_pivot++ )
                    {
                        result += G_tau(i_pivot-N_C-1,j_pivot-N_C) * alpha[m[j_pivot]-1];
                    }
                }
            }

            else
            {
                for ( j_pivot = N_Z ; j_pivot < N ; j_pivot++ )
                {
                    result += get_G_tau(i,m[j_pivot]) * alpha[m[j_pivot]-1];
                }
            }
        }

        if ( !SVM_HINT_RHO_ZERO(svflags) )
        {
            if ( tau[i-1] > 0 )
            {
                result += rho[i-1];
            }

            else if ( tau[i-1] < 0 )
            {
                result -= rho_star[i-1];
            }
        }
    }

    FN_EXIT_POINT result;
}


void SVdata::step_none(L_DOUBLE d_b)
{
    FN_ENTRY_POINT;

    long j_pivot,j_pivot_start;

    if ( E_CACHE_FULL(svflags) )
    {
        j_pivot_start = 1;
    }

    else if ( E_CACHE_SUPP(svflags) )
    {
        j_pivot_start = N_Z+1;
    }

    else if ( E_CACHE_FREE(svflags) )
    {
        j_pivot_start = N_C+1;
    }

    else
    {
        j_pivot_start = N+1;
    }

    if ( j_pivot_start <= N )
    {
        for ( j_pivot = j_pivot_start ; j_pivot <= N ; j_pivot++ )
        {
            e_tau[m[j_pivot-1]-1] += d_b;
        }

        fix_e_free_pivot();
    }

    b += d_b;

    FN_EXIT_POINT;
}

void SVdata::step_one(L_DOUBLE d_alpha, long i)
{
    FN_ENTRY_POINT;

    /*
       NB: we quite deliberately allow for steps when tau = +-2 so-as to
           allow boundary movement
    */

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );
    THROW_ASSERT( tau[i-1] != 0 );

    long j_pivot,j_pivot_start;

    alpha[i-1] += d_alpha;

    if ( !SVM_HINT_MU_ZERO(svflags) )
    {
        if ( tau[i-1] > 0 )
        {
            E += d_alpha * mu[i-1];
        }

        else
        {
            E -= d_alpha * mu_star[i-1];
        }
    }

    f += d_alpha;

    if ( E_CACHE_FULL(svflags) )
    {
        j_pivot_start = 1;
    }

    else if ( E_CACHE_SUPP(svflags) )
    {
        j_pivot_start = N_Z+1;
    }

    else if ( E_CACHE_FREE(svflags) )
    {
        j_pivot_start = N_C+1;
    }

    else
    {
        j_pivot_start = N+1;
    }

    if ( j_pivot_start <= N )
    {
        if ( KERN_CACHE_FULL(svflags) )
        {
            for ( j_pivot = j_pivot_start ; j_pivot <= N ; j_pivot++ )
            {
                e_tau[m[j_pivot-1]-1] += G_tau(i-1,m[j_pivot-1]-1)*d_alpha;
            }
        }

        else if ( KERN_CACHE_FREE(svflags) )
        {
            /*
               Note j_pivot_start = 1,N_Z+1,N_C+1
            */

            long i_pivot;

            i_pivot = minv[i-1];

            /*
               Take care - internal functions can call this for variables
               not technically "free"
            */

            if ( i_pivot <= N_C )
            {
                for ( j_pivot = j_pivot_start ; j_pivot <= N ; j_pivot++ )
                {
                    e_tau[m[j_pivot-1]-1] += get_G_tau(i,m[j_pivot-1])*d_alpha;
                }
            }

            else
            {
                if ( j_pivot_start <= N_C )
                {
                    for ( j_pivot = j_pivot_start ; j_pivot <= N_C ; j_pivot++ )
                    {
                        e_tau[m[j_pivot-1]-1] += get_G_tau(i,m[j_pivot-1])*d_alpha;
                    }
                }

                for ( j_pivot = N_C+1 ; j_pivot <= N ; j_pivot++ )
                {
                    e_tau[m[j_pivot-1]-1] += G_tau(i_pivot-N_C-1,j_pivot-N_C-1)*d_alpha;
                }
            }
        }

        else
        {
            for ( j_pivot = j_pivot_start ; j_pivot <= N ; j_pivot++ )
            {
                e_tau[m[j_pivot-1]-1] += get_G_tau(i,m[j_pivot-1])*d_alpha;
            }
        }

        fix_e_free_pivot();
    }

    FN_EXIT_POINT;
}

void SVdata::step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( d_alpha_pivot.get_effective_size() == d_e_pivot.get_effective_size() );
    THROW_ASSERT( d_alpha_pivot.get_effective_size() == finish_pivot-start_pivot+1 );
    THROW_ASSERT( start_pivot >= N_C+1 );
    THROW_ASSERT( start_pivot <= N );
    THROW_ASSERT( finish_pivot >= N_C+1 );
    THROW_ASSERT( finish_pivot <= N );
    THROW_ASSERT( start_pivot <= finish_pivot );

    long i_pivot,j_pivot,i_pivot_start;

    for ( i_pivot = start_pivot ; i_pivot <= finish_pivot ; i_pivot++ )
    {
        alpha[m[i_pivot-1]-1] += d_alpha_pivot[i_pivot-start_pivot];
    }

    if ( E_CACHE_FULL(svflags) || E_CACHE_SUPP(svflags) || E_CACHE_FREE(svflags) )
    {
        for ( i_pivot = start_pivot+n_ignore ; i_pivot <= finish_pivot ; i_pivot++ )
        {
            e_tau[m[i_pivot-1]-1] += d_e_pivot[i_pivot-start_pivot];
        }
    }

    /*
       update E and f

       E = mu_(star)^T |alpha|
       f = 1^T alpha
    */

    f += d_f;

    if ( !SVM_HINT_MU_ZERO(svflags) )
    {
        for ( i_pivot = start_pivot ; i_pivot <= finish_pivot ; i_pivot++ )
        {
            if ( tau[m[i_pivot-1]-1] > 0 )
            {
                E += d_alpha_pivot[i_pivot-start_pivot] * mu[m[i_pivot-1]-1];
            }

            else
            {
                E -= d_alpha_pivot[i_pivot-start_pivot] * mu_star[m[i_pivot-1]-1];
            }
        }
    }

    /*
       update e_tau
    */

    if ( E_CACHE_FULL(svflags) )
    {
        i_pivot_start = 1;
    }

    else if ( E_CACHE_SUPP(svflags) )
    {
        i_pivot_start = N_Z+1;
    }

    else if ( E_CACHE_FREE(svflags) )
    {
        i_pivot_start = N_C+1;
    }

    else
    {
        i_pivot_start = N+1;
    }

    if ( ( start_pivot > 1 ) && ( i_pivot_start < start_pivot ) )
    {
        for ( i_pivot = i_pivot_start ; i_pivot < start_pivot ; i_pivot++ )
        {
            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    e_tau[m[i_pivot-1]-1] += G_tau(m[i_pivot-1]-1,m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    e_tau[m[i_pivot-1]-1] += get_G_tau(m[i_pivot-1],m[j_pivot-1])*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else if ( KERN_CACHE_DYNA(svflags) )
            {
                /*
                   This cannot happen with the current set of optimisation
                   schemes.
                */

                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    e_tau[m[i_pivot-1]-1] += get_G_tau(m[i_pivot-1],m[j_pivot-1])*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else
            {
                /*
                   This cannot happen with the current set of optimisation
                   schemes.
                */

                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    e_tau[m[i_pivot-1]-1] += get_G_tau(m[i_pivot-1],m[j_pivot-1])*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            e_tau[m[i_pivot-1]-1] += d_b;
        }
    }

    if ( i_pivot_start < finish_pivot+1 )
    {
        i_pivot_start = finish_pivot+1;
    }

    if ( i_pivot_start <= N )
    {
        for ( i_pivot = i_pivot_start ; i_pivot <= N ; i_pivot++ )
        {
            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    e_tau[m[i_pivot-1]-1] += G_tau(m[i_pivot-1]-1,m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    e_tau[m[i_pivot-1]-1] += get_G_tau(m[i_pivot-1],m[j_pivot-1])*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else if ( KERN_CACHE_DYNA(svflags) )
            {
                /*
                   This cannot happen with the current set of optimisation
                   schemes.
                */

                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    e_tau[m[i_pivot-1]-1] += get_G_tau(m[i_pivot-1],m[j_pivot-1])*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else
            {
                /*
                   This cannot happen with the current set of optimisation
                   schemes.
                */

                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    e_tau[m[i_pivot-1]-1] += get_G_tau(m[i_pivot-1],m[j_pivot-1])*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            e_tau[m[i_pivot-1]-1] += d_b;
        }
    }

    b += d_b;

    fix_e_free_pivot();

    FN_EXIT_POINT;
}

void SVdata::step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( d_alpha_pivot.get_effective_size() == finish_pivot-start_pivot+1 );
    THROW_ASSERT( start_pivot >= N_C+1 );
    THROW_ASSERT( start_pivot <= N );
    THROW_ASSERT( finish_pivot >= N_C+1 );
    THROW_ASSERT( finish_pivot <= N );
    THROW_ASSERT( start_pivot <= finish_pivot );

    long i_pivot,j_pivot,i_pivot_start;

    for ( i_pivot = start_pivot ; i_pivot <= finish_pivot ; i_pivot++ )
    {
        alpha[m[i_pivot-1]-1] += d_alpha_pivot[i_pivot-start_pivot];
    }

    /*
       update E and f

       E = mu_(star)^T |alpha|
       f = 1^T alpha
    */

    if ( !SVM_HINT_MU_ZERO(svflags) )
    {
        for ( i_pivot = start_pivot ; i_pivot <= finish_pivot ; i_pivot++ )
        {
            if ( tau[m[i_pivot-1]-1] > 0 )
            {
                E += d_alpha_pivot[i_pivot-start_pivot] * mu[m[i_pivot-1]-1];
            }

            else
            {
                E -= d_alpha_pivot[i_pivot-start_pivot] * mu_star[m[i_pivot-1]-1];
            }
        }
    }

    for ( i_pivot = start_pivot ; i_pivot <= finish_pivot ; i_pivot++ )
    {
        f += d_alpha_pivot[i_pivot-start_pivot];
    }

    /*
       update e_tau
    */

    if ( E_CACHE_FULL(svflags) )
    {
        i_pivot_start = 1;
    }

    else if ( E_CACHE_SUPP(svflags) )
    {
        i_pivot_start = N_Z+1;
    }

    else if ( E_CACHE_FREE(svflags) )
    {
        i_pivot_start = N_C+1;
    }

    else
    {
        i_pivot_start = N+1;
    }

    if ( i_pivot_start <= N )
    {
        for ( i_pivot = i_pivot_start ; i_pivot <= N ; i_pivot++ )
        {
            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    e_tau[m[i_pivot-1]-1] += G_tau(m[i_pivot-1]-1,m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    e_tau[m[i_pivot-1]-1] += get_G_tau(m[i_pivot-1],m[j_pivot-1])*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else if ( KERN_CACHE_DYNA(svflags) )
            {
                /*
                   This cannot happen with the current set of optimisation
                   schemes.
                */

                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    e_tau[m[i_pivot-1]-1] += get_G_tau(m[i_pivot-1],m[j_pivot-1])*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else
            {
                /*
                   This cannot happen with the current set of optimisation
                   schemes.
                */

                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    e_tau[m[i_pivot-1]-1] += get_G_tau(m[i_pivot-1],m[j_pivot-1])*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            e_tau[m[i_pivot-1]-1] += d_b;
        }

        fix_e_free_pivot();
    }

    b += d_b;

    FN_EXIT_POINT;
}

L_DOUBLE SVdata::get_G_no_extras(long i, long j)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );
    THROW_ASSERT( j >= 1 );
    THROW_ASSERT( j <= N );

    L_DOUBLE result = 0.0;

    if ( KERN_CACHE_FULL(svflags) || KERN_CACHE_FREE(svflags) )
    {
        result = get_G_tau(i,j);

        tryagain:

        if ( !SVM_HINT_GAMMA_ZERO(svflags) )
        {
            if ( i == j )
            {
                if ( tau[i-1] > 0 )
                {
                    result -= gamma[i-1];
                }

                else if ( tau[i-1] < 0 )
                {
                    result -= gamma_star[i-1];
                }
            }
        }

        if ( !SVM_HINT_MU_ZERO(svflags) )
        {
            if ( ( tau[i-1] > 0 ) && ( tau[j-1] > 0 ) )
            {
                result -= mu[i-1] * mu[j-1];
            }

            else if ( ( tau[i-1] > 0 ) && ( tau[j-1] < 0 ) )
            {
                result += mu[i-1] * mu_star[j-1];
            }

            else if ( ( tau[i-1] < 0 ) && ( tau[j-1] > 0 ) )
            {
                result += mu_star[i-1] * mu[j-1];
            }

            else if ( ( tau[i-1] < 0 ) && ( tau[j-1] < 0 ) )
            {
                result -= mu_star[i-1] * mu_star[j-1];
            }
        }
    }

    else if ( KERN_CACHE_DYNA(svflags) )
    {
        result = kerncache.getval_statcache(this,i,j);
    }

    else
    {
        result = calc_G_tau_from_scratch(i,j);

        goto tryagain;
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::get_G_tau(long i, long j)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );
    THROW_ASSERT( j >= 1 );
    THROW_ASSERT( j <= N );

    L_DOUBLE result;
    long i_pivot,j_pivot;

    if ( KERN_CACHE_FULL(svflags) )
    {
        result = G_tau(i-1,j-1);
    }

    else if ( KERN_CACHE_FREE(svflags) )
    {
        i_pivot = minv[i-1];
        j_pivot = minv[j-1];

        if ( ( i_pivot > N_C ) && ( j_pivot > N_C ) )
        {
            result = G_tau(i_pivot-N_C-1,j_pivot-N_C-1);
        }

        else
        {
            result = calc_G_tau_from_scratch(i,j);
        }
    }

    else if ( KERN_CACHE_DYNA(svflags) )
    {
        result = kerncache.getval_statcache(this,i,j);

        if ( !SVM_HINT_GAMMA_ZERO(svflags) )
        {
            if ( i == j )
            {
                if ( tau[i-1] > 0 )
                {
                    result += gamma[i-1];
                }

                else if ( tau[i-1] < 0 )
                {
                    result += gamma_star[i-1];
                }
            }
        }

        if ( !SVM_HINT_MU_ZERO(svflags) )
        {
            if ( ( tau[i-1] > 0 ) && ( tau[j-1] > 0 ) )
            {
                result += mu[i-1] * mu[j-1];
            }

            else if ( ( tau[i-1] > 0 ) && ( tau[j-1] < 0 ) )
            {
                result -= mu[i-1] * mu_star[j-1];
            }

            else if ( ( tau[i-1] < 0 ) && ( tau[j-1] > 0 ) )
            {
                result -= mu_star[i-1] * mu[j-1];
            }

            else if ( ( tau[i-1] < 0 ) && ( tau[j-1] < 0 ) )
            {
                result += mu_star[i-1] * mu_star[j-1];
            }
        }
    }

    else
    {
        result = calc_G_tau_from_scratch(i,j);
    }

    FN_EXIT_POINT result;
}

void SVdata::constrain_Z(long i)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );

    constrain_Z_pivot(minv[i-1]);

    FN_EXIT_POINT;
}

void SVdata::constrain_L(long i)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );

    constrain_L_pivot(minv[i-1]);

    FN_EXIT_POINT;
}

void SVdata::constrain_U(long i)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );

    constrain_U_pivot(minv[i-1]);

    FN_EXIT_POINT;
}

void SVdata::free_L(long i)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );

    free_L_pivot(minv[i-1]);

    FN_EXIT_POINT;
}

void SVdata::free_U(long i)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );

    free_U_pivot(minv[i-1]);

    FN_EXIT_POINT;
}

void SVdata::constrain_Z_pivot(long i_pivot)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivot >= N_C+1 );
    THROW_ASSERT( i_pivot <= N );
    THROW_ASSERT( tau[m[i_pivot-1]-1] != -2 );
    THROW_ASSERT( tau[m[i_pivot-1]-1] != +2 );

    long i,j,j_pivot;

    if ( tau[m[i_pivot-1]-1] != 0 )
    {
        /*
           fix e_free_pivot
        */

        e_free_pivot.remove(i_pivot-N_C);

        i = m[i_pivot-1];

        if ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) )
        {
            /*
               Fixup beth_tau
            */

            SHRINK_BETH(i_pivot-N_C);
        }

        if ( KERN_CACHE_FREE(svflags) )
        {
            G_tau.remove(i_pivot-N_C);
        }

        if ( tau[i-1] == -1 )
        {
            /*
               Fixup G_tau - this operation takes advantage of the symmetry
               property of the hessian matrix class, so that the operation
               must be done for only one of the row/column.
            */

            if ( KERN_CACHE_FULL(svflags) )
            {
                if ( !SVM_HINT_GAMMA_ZERO(svflags) )
                {
                    for ( j_pivot = N_Z+1 ; j_pivot <= N ; j_pivot++ )
                    {
                        j = m[j_pivot-1];

                        if ( i == j )
                        {
                            G_tau(i-1,j-1) -= gamma_star[i-1];
                        }
                    }
                }

                if ( !SVM_HINT_MU_ZERO(svflags) )
                {
                    for ( j_pivot = N_Z+1 ; j_pivot <= N ; j_pivot++ )
                    {
                        j = m[j_pivot-1];

                        if ( tau[j-1] < 0 )
                        {
                            G_tau(i-1,j-1) -= ( mu_star[i-1] * mu_star[j-1] );
                        }

                        else
                        {
                            G_tau(i-1,j-1) += ( mu_star[i-1] * mu[j-1] );
                        }
                    }
                }
            }

            /*
               Fixup e_tau
            */

            if ( E_CACHE_FULL(svflags) && !SVM_HINT_RHO_ZERO(svflags) && ( contype[i-1] != 1 ) )
            {
                e_tau[i-1] += rho_star[i-1];
            }

            if ( E_CACHE_FULL(svflags) && !SVM_HINT_MU_ZERO(svflags) )
            {
                e_tau[i-1] += E * mu_star[i-1];
            }

            N_FN--;
        }

        else
        {
            /*
               Fixup G_tau - this operation takes advantage of the symmetry
               property of the hessian matrix class, so that the operation
               must be done for only one of the row/column.
            */

            if ( KERN_CACHE_FULL(svflags) )
            {
                if ( !SVM_HINT_GAMMA_ZERO(svflags) )
                {
                    for ( j_pivot = N_Z+1 ; j_pivot <= N ; j_pivot++ )
                    {
                        j = m[j_pivot-1];

                        if ( i == j )
                        {
                            G_tau(i-1,j-1) -= gamma[i-1];
                        }
                    }
                }

                if ( !SVM_HINT_MU_ZERO(svflags) )
                {
                    for ( j_pivot = N_Z+1 ; j_pivot <= N ; j_pivot++ )
                    {
                        j = m[j_pivot-1];

                        if ( tau[j-1] < 0 )
                        {
                            G_tau(i-1,j-1) += ( mu[i-1] * mu_star[j-1] );
                        }

                        else
                        {
                            G_tau(i-1,j-1) -= ( mu[i-1] * mu[j-1] );
                        }
                    }
                }
            }

            /*
               Fixup e_tau
            */

            if ( E_CACHE_FULL(svflags) && !SVM_HINT_RHO_ZERO(svflags) && ( contype[i-1] != 2 ) )
            {
                e_tau[i-1] -= rho[i-1];
            }

            if ( E_CACHE_FULL(svflags) && !SVM_HINT_MU_ZERO(svflags) )
            {
                e_tau[i-1] -= E * mu[i-1];
            }

            N_FP--;
        }

        /*
           Fixup tau
        */

        tau[i-1] = 0;

        /*
           Fixup the order vector
        */

        m.bswap(N_Z+1,i_pivot);

        for ( j_pivot = N_Z+1 ; j_pivot <= i_pivot ; j_pivot++ )
        {
            minv[m[j_pivot-1]-1] = j_pivot;
        }

        N_Z++;
        N_F--;
        N_C++;
        N_S--;
    }

    FN_EXIT_POINT;
}

void SVdata::constrain_L_pivot(long i_pivot)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivot >= N_C+1 );
    THROW_ASSERT( i_pivot <= N );
    THROW_ASSERT( tau[m[i_pivot-1]-1] != 0 );
    THROW_ASSERT( tau[m[i_pivot-1]-1] != +1 );
    THROW_ASSERT( tau[m[i_pivot-1]-1] != +2 );

    long i,j_pivot;

    if ( tau[m[i_pivot-1]-1] != -2 )
    {
        /*
           fix e_free_pivot
        */

        e_free_pivot.remove(i_pivot-N_C);

        i = m[i_pivot-1];

        if ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) )
        {
            /*
               Fixup beth_tau
            */

            SHRINK_BETH(i_pivot-N_C);
        }

        if ( KERN_CACHE_FREE(svflags) )
        {
            G_tau.remove(i_pivot-N_C);
        }

        /*
           Neither G_tau nor e_tau will be changed by this operation.
        */

        /*
           Fixup tau
        */

        tau[i-1] = -2;

        /*
           Fixup the order vector
        */

        m.bswap(N_Z+N_L+1,i_pivot);

        for ( j_pivot = N_Z+N_L+1 ; j_pivot <= i_pivot ; j_pivot++ )
        {
            minv[m[j_pivot-1]-1] = j_pivot;
        }

        N_L++;
        N_F--;
        N_FN--;
        N_C++;
    }

    FN_EXIT_POINT;
}

void SVdata::constrain_U_pivot(long i_pivot)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivot >= N_C+1 );
    THROW_ASSERT( i_pivot <= N );
    THROW_ASSERT( tau[m[i_pivot-1]-1] != -2 );
    THROW_ASSERT( tau[m[i_pivot-1]-1] != -1 );
    THROW_ASSERT( tau[m[i_pivot-1]-1] != 0 );

    long i,j_pivot;

    if ( tau[m[i_pivot-1]-1] != +2 )
    {
        /*
           fix e_free_pivot
        */

        e_free_pivot.remove(i_pivot-N_C);

        i = m[i_pivot-1];

        if ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) )
        {
            /*
               Fixup beth_tau
            */

            SHRINK_BETH(i_pivot-N_C);
        }

        if ( KERN_CACHE_FREE(svflags) )
        {
            G_tau.remove(i_pivot-N_C);
        }

        /*
           Neither G_tau nor e_tau will be changed by this operation.
        */

        /*
           Fixup tau
        */

        tau[i-1] = +2;

        /*
           Fixup the order vector
        */

        m.bswap(N_C+1,i_pivot);

        for ( j_pivot = N_C+1 ; j_pivot <= i_pivot ; j_pivot++ )
        {
            minv[m[j_pivot-1]-1] = j_pivot;
        }

        N_U++;
        N_F--;
        N_FP--;
        N_C++;
    }

    FN_EXIT_POINT;
}

void SVdata::free_L_pivot(long i_pivot)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivot >= 1 );
    THROW_ASSERT( i_pivot <= N_Z+N_L+1 );
    THROW_ASSERT( tau[m[i_pivot-1]-1] != 1 );
    THROW_ASSERT( tau[m[i_pivot-1]-1] != 2 );
    THROW_ASSERT( contype[m[i_pivot-1]-1] != 0 );
    THROW_ASSERT( contype[m[i_pivot-1]-1] != 2 );

    long i,j,j_pivot;
    int old_tau;

    if ( tau[m[i_pivot-1]-1] != -1 )
    {
        i = m[i_pivot-1];

        old_tau = tau[i-1];

        /*
           If e is entering the cache, will need to ensure that it has a
           valid value before we start.
        */

        if ( ( E_CACHE_SUPP(svflags) && ( old_tau == 0 ) ) || ( E_CACHE_FREE(svflags) ) )
        {
            e_tau[i-1] = get_e(i);

            if ( contype[i-1] == 1 )
            {
                e_tau[i-1] -= rho_star[i-1];
            }
        }

        /*
           Pivotting will be needed first.
        */

        if ( old_tau == 0 )
        {
            N_Z--;
            N_F++;
            N_FN++;
            N_C--;
            N_S++;
        }

        else
        {
            N_L--;
            N_F++;
            N_FN++;
            N_C--;
        }

        /*
           Fixup tau
        */

        tau[i-1] = -1;

        /*
           Fixup the order vector
        */

        m.fswap(i_pivot,N);

        for ( j_pivot = i_pivot ; j_pivot <= N ; j_pivot++ )
        {
            minv[m[j_pivot-1]-1] = j_pivot;
        }

        /*
           At this point, the point has moved, so reset i_pivot.
        */

        i_pivot = N;

        if ( old_tau == 0 )
        {
            /*
               Fixup G_tau - this operation takes advantage of the symmetry
               property of the hessian matrix class, so that the operation
               must be done for only one of the row/column.
            */

            if ( KERN_CACHE_FULL(svflags) )
            {
                if ( !SVM_HINT_GAMMA_ZERO(svflags) )
                {
                    for ( j_pivot = N_Z+1 ; j_pivot <= N ; j_pivot++ )
                    {
                        j = m[j_pivot-1];

                        if ( i == j )
                        {
                            G_tau(i-1,j-1) += gamma_star[i-1];
                        }
                    }
                }

                if ( !SVM_HINT_MU_ZERO(svflags) )
                {
                    for ( j_pivot = N_Z+1 ; j_pivot <= N ; j_pivot++ )
                    {
                        j = m[j_pivot-1];

                        if ( tau[j-1] < 0 )
                        {
                            G_tau(i-1,j-1) += ( mu_star[i-1] * mu_star[j-1] );
                        }

                        else
                        {
                            G_tau(i-1,j-1) -= ( mu_star[i-1] * mu[j-1] );
                        }
                    }
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                fVECTOR ax('x',N_F);

                for ( j_pivot = N_C+1 ; j_pivot <= N ; j_pivot++ )
                {
                    ax[j_pivot-N_C-1] = calc_G_tau_from_scratch(m[i_pivot-1],m[j_pivot-1]);
                }

                G_tau.addend(ax[N_F-1],ax);
            }

            /*
               Fixup e_tau
            */

            if ( ( E_CACHE_FULL(svflags) || E_CACHE_SUPP(svflags) || E_CACHE_FREE(svflags) ) && !SVM_HINT_RHO_ZERO(svflags) && ( ( contype[i-1] != 1 ) || ( contype[i-1] != 1 ) ) )
            {
                e_tau[i-1] -= rho_star[i-1];
            }

            if ( ( E_CACHE_FULL(svflags) || E_CACHE_SUPP(svflags) || E_CACHE_FREE(svflags) ) && !SVM_HINT_MU_ZERO(svflags) )
            {
                e_tau[i-1] -= E * mu_star[i-1];
            }
        }

        else if ( KERN_CACHE_FREE(svflags) )
        {
            fVECTOR ax('x',N_F);

            for ( j_pivot = N_C+1 ; j_pivot <= N ; j_pivot++ )
            {
                ax[j_pivot-N_C-1] = calc_G_tau_from_scratch(m[i_pivot-1],m[j_pivot-1]);
            }

            G_tau.addend(ax[N_F-1],ax);
        }

        /*
           Add new element to e_free_pivot
        */

        e_free_pivot.addend(get_e(i));

        if ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) )
        {
            /*
               Fixup beth_tau
            */

            fVECTOR ax('x',N_F);

            /*
               Implicitly, if a factorisation is kept then G_tau must
               be either FULL or FREE
            */

            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j_pivot = N_C+1 ; j_pivot <= N ; j_pivot++ )
                {
                    ax[j_pivot-N_C-1] = G_tau(m[i_pivot-1]-1,m[j_pivot-1]-1);
                }
            }

            else
            {
                for ( j_pivot = N_C+1 ; j_pivot <= N ; j_pivot++ )
                {
                    ax[j_pivot-N_C-1] = G_tau(i_pivot-N_C-1,j_pivot-N_C-1);
                }
            }

            ADDEND_BETH(ax);
        }
    }

    FN_EXIT_POINT;
}

void SVdata::free_U_pivot(long i_pivot)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( ( ( i_pivot >= 1 ) && ( i_pivot <= N_Z ) ) || ( ( i_pivot >= N_Z+N_L+1 ) && ( i_pivot <= N_C ) ) );
    THROW_ASSERT( tau[m[i_pivot-1]-1] != -1 );
    THROW_ASSERT( tau[m[i_pivot-1]-1] != -2 );
    THROW_ASSERT( contype[m[i_pivot-1]-1] != 0 );
    THROW_ASSERT( contype[m[i_pivot-1]-1] != 1 );

    long i,j,j_pivot;
    int old_tau;

    if ( tau[m[i_pivot-1]-1] != +1 )
    {
        i = m[i_pivot-1];

        old_tau = tau[i-1];

        /*
           If e is entering the cache, will need to ensure that it has a
           valid value before we start.
        */

        if ( ( E_CACHE_SUPP(svflags) && ( old_tau == 0 ) ) || ( E_CACHE_FREE(svflags) ) )
        {
            e_tau[i-1] = get_e(i);

            if ( contype[i-1] == 2 )
            {
                e_tau[i-1] += rho[i-1];
            }
        }

        if ( old_tau == 0 )
        {
            N_Z--;
            N_F++;
            N_FP++;
            N_C--;
            N_S++;
        }

        else
        {
            N_U--;
            N_F++;
            N_FP++;
            N_C--;
        }

        /*
           Fixup tau
        */

        tau[i-1] = +1;

        /*
           Fixup the order vector
        */

        m.fswap(i_pivot,N);

        for ( j_pivot = i_pivot ; j_pivot <= N ; j_pivot++ )
        {
            minv[m[j_pivot-1]-1] = j_pivot;
        }

        /*
           At this point, the point has moved, so reset i_pivot.
        */

        i_pivot = N;

        if ( old_tau == 0 )
        {
            /*
               Fixup G_tau - this operation takes advantage of the symmetry
               property of the hessian matrix class, so that the operation
               must be done for only one of the row/column.
            */

            if ( KERN_CACHE_FULL(svflags) )
            {
                if ( !SVM_HINT_GAMMA_ZERO(svflags) )
                {
                    for ( j_pivot = N_Z+1 ; j_pivot <= N ; j_pivot++ )
                    {
                        j = m[j_pivot-1];

                        if ( i == j )
                        {
                            G_tau(i-1,j-1) += gamma[i-1];;
                        }
                    }
                }

                if ( !SVM_HINT_MU_ZERO(svflags) )
                {
                    for ( j_pivot = N_Z+1 ; j_pivot <= N ; j_pivot++ )
                    {
                        j = m[j_pivot-1];

                        if ( tau[j-1] < 0 )
                        {
                            G_tau(i-1,j-1) -= ( mu[i-1] * mu_star[j-1] );
                        }

                        else
                        {
                            G_tau(i-1,j-1) += ( mu[i-1] * mu[j-1] );
                        }
                    }
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                fVECTOR ax('x',N_F);

                for ( j_pivot = N_C+1 ; j_pivot <= N ; j_pivot++ )
                {
                    ax[j_pivot-N_C-1] = calc_G_tau_from_scratch(m[i_pivot-1],m[j_pivot-1]);
                }

                G_tau.addend(ax[N_F-1],ax);
            }

            /*
               Fixup e_tau
            */

            if ( ( E_CACHE_FULL(svflags) || E_CACHE_SUPP(svflags) || E_CACHE_FREE(svflags) ) && !SVM_HINT_RHO_ZERO(svflags) && ( contype[i-1] != 2 ) )
            {
                e_tau[i-1] += rho[i-1];
            }

            if ( ( E_CACHE_FULL(svflags) || E_CACHE_SUPP(svflags) || E_CACHE_FREE(svflags) ) && !SVM_HINT_MU_ZERO(svflags) )
            {
                e_tau[i-1] += E * mu[i-1];
            }
        }

        else if ( KERN_CACHE_FREE(svflags) )
        {
            fVECTOR ax('x',N_F);

            for ( j_pivot = N_C+1 ; j_pivot <= N ; j_pivot++ )
            {
                ax[j_pivot-N_C-1] = calc_G_tau_from_scratch(m[i_pivot-1],m[j_pivot-1]);
            }

            G_tau.addend(ax[N_F-1],ax);
        }

        /*
           Add new element to e_free_pivot
        */

        e_free_pivot.addend(get_e(i));

        if ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) )
        {
            /*
               Fixup beth_tau
            */

            fVECTOR ax('x',N_F);

            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j_pivot = N_C+1 ; j_pivot <= N ; j_pivot++ )
                {
                    ax[j_pivot-N_C-1] = G_tau(m[i_pivot-1]-1,m[j_pivot-1]-1);
                }
            }

            else
            {
                for ( j_pivot = N_C+1 ; j_pivot <= N ; j_pivot++ )
                {
                    ax[j_pivot-N_C-1] = G_tau(i_pivot-N_C-1,j_pivot-N_C-1);
                }
            }

            ADDEND_BETH(ax);
        }
    }

    FN_EXIT_POINT;
}

long SVdata::find_marker(long id)
{
    FN_ENTRY_POINT;

    long result = 0;
    long i;

    if ( N > 0 )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            if ( ident[i-1] == id )
            {
                result = i;
            }
        }
    }

    FN_EXIT_POINT result;
}

int SVdata::is_sparse(void) const
{
    FN_ENTRY_POINT;

    int result;

    result = K.is_sparse();

    FN_EXIT_POINT result;
}

void SVdata::fix_e_free_pivot(void)
{
    FN_ENTRY_POINT;

    long j_pivot;

    if ( N_F >= 1 )
    {
        for ( j_pivot = N_C+1 ; j_pivot <= N ; j_pivot++ )
        {
            e_free_pivot[j_pivot-N_C-1] = get_e(m[j_pivot-1]);
        }
    }

    FN_EXIT_POINT;
}


void SVdata::regen_from_short_summary(const SVdata_short_summary &source)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( N == source.N );

    e_free_pivot.make_zero(DO_FORCE);
    e_free_pivot.make_normal(DO_FORCE);

    contype      = source.contype;
    alpha        = source.alpha;
    tau          = source.tau;
    m            = source.m;
    minv         = source.minv;
    e_tau        = source.e_tau;
    e_free_pivot = source.e_free_pivot;

    b = source.b;
    E = source.E;
    f = source.f;

    beth_tau = source.beth_tau;
    svflags  = source.svflags;

    N    = source.N;
    N_Z  = source.N_Z;
    N_L  = source.N_L;
    N_U  = source.N_U;
    N_F  = source.N_F;
    N_FN = source.N_FN;
    N_FP = source.N_FP;
    N_C  = source.N_C;
    N_S  = source.N_S;

    cache_memsize    = source.cache_memsize;
    cache_min_rowdim = source.cache_min_rowdim;

    d2c_buff_size = source.d2c_buff_size;

    if ( KERN_CACHE_DYNA(svflags) )
    {
        kerncache.reset();
        kerncache.setmemsize(this,cache_memsize,cache_min_rowdim);
    }

    FN_EXIT_POINT;
}

SVdata_short_summary::SVdata_short_summary(const SVdata &source)
{
    FN_ENTRY_POINT;

    contype.make_zero(DO_FORCE);
    alpha.make_zero(DO_FORCE);
    tau.make_zero(DO_FORCE);
    m.make_zero(DO_FORCE);
    minv.make_zero(DO_FORCE);
    e_tau.make_zero(DO_FORCE);
    e_free_pivot.make_zero(DO_FORCE);

    contype      = source.contype;
    alpha        = source.alpha;
    tau          = source.tau;
    m            = source.m;
    minv         = source.minv;
    e_tau        = source.e_tau;
    e_free_pivot = source.e_free_pivot;

    b = source.b;
    E = source.E;
    f = source.f;

    beth_tau = source.beth_tau;
    svflags  = source.svflags;

    N    = source.N;
    N_Z  = source.N_Z;
    N_L  = source.N_L;
    N_U  = source.N_U;
    N_F  = source.N_F;
    N_FN = source.N_FN;
    N_FP = source.N_FP;
    N_C  = source.N_C;
    N_S  = source.N_S;

    cache_memsize    = source.cache_memsize;
    cache_min_rowdim = source.cache_min_rowdim;

    d2c_buff_size = source.d2c_buff_size;

    FN_EXIT_POINT;
}

SVdata_summary::SVdata_summary(const SVdata &source)
{
    FN_ENTRY_POINT;

    long i;

    fix_bias = source.fix_bias;
    svflags  = source.svflags;

    N    = source.N;
    N_Z  = source.N_Z;
    N_L  = source.N_L;
    N_U  = source.N_U;
    N_F  = source.N_F;
    N_FP = source.N_FP;
    N_FN = source.N_FN;
    N_C  = source.N_C;
    N_S  = source.N_S;

    b = source.b;
    E = source.E;
    f = source.f;

    x.trim_to_size(0);

    if ( N > 0 )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            x.addend((source.x).get_offset_element(i-1));
        }
    }

    z.make_zero(DO_FORCE);             z            = source.z;
    rho.make_zero(DO_FORCE);           rho          = source.rho;
    rho_star.make_zero(DO_FORCE);      rho_star     = source.rho_star;
    v.make_zero(DO_FORCE);             v            = source.v;
    h.make_zero(DO_FORCE);             h            = source.h;
    gamma.make_zero(DO_FORCE);         gamma        = source.gamma;
    gamma_star.make_zero(DO_FORCE);    gamma_star   = source.gamma_star;
    mu.make_zero(DO_FORCE);            mu           = source.mu;
    mu_star.make_zero(DO_FORCE);       mu_star      = source.mu_star;
    contype.make_zero(DO_FORCE);       contype      = source.contype;
    ident.make_zero(DO_FORCE);         ident        = source.ident;
    free_ones.make_zero(DO_FORCE);     free_ones    = source.free_ones;
    alpha.make_zero(DO_FORCE);         alpha        = source.alpha;
    tau.make_zero(DO_FORCE);           tau          = source.tau;
    m.make_zero(DO_FORCE);             m            = source.m;
    minv.make_zero(DO_FORCE);          minv         = source.minv;
    e_tau.make_zero(DO_FORCE);         e_tau        = source.e_tau;
                                       beth_tau     = source.beth_tau;
    e_free_pivot.make_zero(DO_FORCE);  e_free_pivot = source.e_free_pivot;

    K = source.K;

    cache_memsize    = source.cache_memsize;
    cache_min_rowdim = source.cache_min_rowdim;

    d2c_buff_size = source.d2c_buff_size;

    FN_EXIT_POINT;
}

void SVdata::remake_from_summary(const SVdata_summary &source)
{
    FN_ENTRY_POINT;

    long i;

    fix_bias = source.fix_bias;
    svflags  = source.svflags;

    N    = source.N;
    N_Z  = source.N_Z;
    N_L  = source.N_L;
    N_U  = source.N_U;
    N_F  = source.N_F;
    N_FP = source.N_FP;
    N_FN = source.N_FN;
    N_C  = source.N_C;
    N_S  = source.N_S;

    b = source.b;
    E = source.E;
    f = source.f;

    x.trim_to_size(0);

    if ( N > 0 )
    {
        for ( i = 1 ; i <= N ; i++ )
        {
            x.addend((source.x).get_offset_element(i-1));
        }
    }

    z.make_zero(DO_FORCE);             z            = source.z;
    rho.make_zero(DO_FORCE);           rho          = source.rho;
    rho_star.make_zero(DO_FORCE);      rho_star     = source.rho_star;
    v.make_zero(DO_FORCE);             v            = source.v;
    h.make_zero(DO_FORCE);             h            = source.h;
    gamma.make_zero(DO_FORCE);         gamma        = source.gamma;
    gamma_star.make_zero(DO_FORCE);    gamma_star   = source.gamma_star;
    mu.make_zero(DO_FORCE);            mu           = source.mu;
    mu_star.make_zero(DO_FORCE);       mu_star      = source.mu_star;
    contype.make_zero(DO_FORCE);       contype      = source.contype;
    ident.make_zero(DO_FORCE);         ident        = source.ident;
    free_ones.make_zero(DO_FORCE);     free_ones    = source.free_ones;
    alpha.make_zero(DO_FORCE);         alpha        = source.alpha;
    tau.make_zero(DO_FORCE);           tau          = source.tau;
    m.make_zero(DO_FORCE);             m            = source.m;
    minv.make_zero(DO_FORCE);          minv         = source.minv;
    e_tau.make_zero(DO_FORCE);         e_tau        = source.e_tau;
                                       beth_tau     = source.beth_tau;
    e_free_pivot.make_zero(DO_FORCE);  e_free_pivot = source.e_free_pivot;

    K = source.K;

    cache_memsize    = source.cache_memsize;
    cache_min_rowdim = source.cache_min_rowdim;

    d2c_buff_size = source.d2c_buff_size;

    if ( KERN_CACHE_DYNA(svflags) )
    {
        kerncache.reset();
        kerncache.setmemsize(this,cache_memsize,cache_min_rowdim);
    }

    /*
       Remake G_tau
    */

    rewrite_G_tau(1);

    FN_EXIT_POINT;
}

L_DOUBLE SVdata::calc_G_tau_from_scratch(long i, long j)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i >= 1 );
    THROW_ASSERT( i <= N );
    THROW_ASSERT( j >= 1 );
    THROW_ASSERT( j <= N );

    L_DOUBLE result;

    result = K.kernel(x[i-1],x[j-1],NULL,NULL,0.0,0.0,i,j,z[i-1],z[j-1]);

    if ( tau[i-1] < 0 )
    {
        if ( !SVM_HINT_GAMMA_ZERO(svflags) )
        {
            if ( i == j )
            {
                result += gamma_star[i-1];
            }
        }

        if ( !SVM_HINT_MU_ZERO(svflags) )
        {
            if ( tau[j-1] < 0 )
            {
                result += mu_star[i-1]*mu_star[j-1];
            }

            else if ( tau[j-1] > 0 )
            {
                result -= mu_star[i-1]*mu[j-1];
            }
        }
    }

    else if ( tau[i-1] > 0 )
    {
        if ( !SVM_HINT_GAMMA_ZERO(svflags) )
        {
            if ( i == j )
            {
                result += gamma[i-1];
            }
        }

        if ( !SVM_HINT_MU_ZERO(svflags) )
        {
            if ( tau[j-1] < 0 )
            {
                result -= mu[i-1]*mu_star[j-1];
            }

            else if ( tau[j-1] > 0 )
            {
                result += mu[i-1]*mu[j-1];
            }
        }
    }

    FN_EXIT_POINT result;
}

void SVdata::rewrite_G_tau(int pad_first)
{
    FN_ENTRY_POINT;

    long i,j,i_pivot,j_pivot;

    /*
       fix G_tau - important to remember matrix symmetry here (don't want
       to add anything twice by mistake).
    */

    if ( KERN_CACHE_FULL(svflags) )
    {
        if ( pad_first )
        {
            G_tau.make_zero(DO_FORCE);
            G_tau.make_symmetric(DO_FORCE);
            G_tau.pad_matrix(N);
        }

        if ( N > 0 )
        {
            for ( i = 1 ; i <= N ; i++ )
            {
                for ( j = 1 ; j <= i ; j++ )
                {
                    G_tau(i-1,j-1) = calc_G_tau_from_scratch(i,j);
                }
            }
        }
    }

    else if ( KERN_CACHE_FREE(svflags) )
    {
        if ( pad_first )
        {
            G_tau.make_zero(DO_FORCE);
            G_tau.make_symmetric(DO_FORCE);
            G_tau.pad_matrix(N_F);
        }

        if ( N_F > 0 )
        {
            for ( i_pivot = N_C+1 ; i_pivot <= N ; i_pivot++ )
            {
                for ( j_pivot = N_C+1 ; j_pivot <= i_pivot ; j_pivot++ )
                {
                    G_tau(i_pivot-N_C-1,j_pivot-N_C-1) = calc_G_tau_from_scratch(minv[i_pivot-1],minv[j_pivot-1]);
                }
            }
        }
    }

    else if ( KERN_CACHE_DYNA(svflags) )
    {
        kerncache.flush();
    }

    FN_EXIT_POINT;
}

std::ostream &operator<<(std::ostream &output, SVdata &dumpee)
{
    FN_ENTRY_POINT;

    output << "SVM data type: N\n\n";

    output << "Bias type: " << dumpee.fix_bias << "\n";
    output << "Flags:     " << dumpee.svflags  << "\n\n";

    output << "N:    " << dumpee.N    << "\n";
    output << "N_Z:  " << dumpee.N_Z  << "\n";
    output << "N_L:  " << dumpee.N_L  << "\n";
    output << "N_U:  " << dumpee.N_U  << "\n";
    output << "N_F:  " << dumpee.N_F  << "\n";
    output << "N_FP: " << dumpee.N_FP << "\n";
    output << "N_FN: " << dumpee.N_FN << "\n";
    output << "N_C:  " << dumpee.N_C  << "\n";
    output << "N_S:  " << dumpee.N_S  << "\n\n";

    output << "b: " << fixnum(dumpee.b) << "\n";
    output << "E: " << fixnum(dumpee.E) << "\n";
    output << "f: " << fixnum(dumpee.f) << "\n\n";

    output << "x:            " << dumpee.x            << "\n";
    output << "z:            " << dumpee.z            << "\n";
    output << "rho:          " << dumpee.rho          << "\n";
    output << "rho_star:     " << dumpee.rho_star     << "\n";
    output << "v:            " << dumpee.v            << "\n";
    output << "h:            " << dumpee.h            << "\n";
    output << "gamma:        " << dumpee.gamma        << "\n";
    output << "gamma_star:   " << dumpee.gamma_star   << "\n";
    output << "mu:           " << dumpee.mu           << "\n";
    output << "mu_star:      " << dumpee.mu_star      << "\n";
    output << "contype:      " << dumpee.contype      << "\n";
    output << "ident:        " << dumpee.ident        << "\n";
    output << "free_ones:    " << dumpee.free_ones    << "\n";
    output << "alpha:        " << dumpee.alpha        << "\n";
    output << "tau:          " << dumpee.tau          << "\n";
    output << "m:            " << dumpee.m            << "\n";
    output << "minv:         " << dumpee.minv         << "\n";
    output << "e_tau:        " << dumpee.e_tau        << "\n";
    output << "e_free_pivot: " << dumpee.e_free_pivot << "\n\n";

    output << "beth_tau: " << dumpee.beth_tau << "\n\n";

    output << "Cache size (MB):     " << dumpee.cache_memsize    << "\n";
    output << "Min cache dimension: " << dumpee.cache_min_rowdim << "\n";
    output << "d2c buffer size:     " << dumpee.d2c_buff_size    << "\n\n";

    output << "Kernel: " << dumpee.K << "\n\n";

    output << "G_tau: " << dumpee.G_tau << "\n\n";

    FN_EXIT_POINT output;
}

std::ostream &operator<<(std::ostream &output, SVdata_summary &dumpee)
{
    FN_ENTRY_POINT;

    output << "SVM data type: S\n\n";

    output << "Bias type: " << dumpee.fix_bias << "\n";
    output << "Flags:     " << dumpee.svflags  << "\n\n";

    output << "N:    " << dumpee.N    << "\n";
    output << "N_Z:  " << dumpee.N_Z  << "\n";
    output << "N_L:  " << dumpee.N_L  << "\n";
    output << "N_U:  " << dumpee.N_U  << "\n";
    output << "N_F:  " << dumpee.N_F  << "\n";
    output << "N_FP: " << dumpee.N_FP << "\n";
    output << "N_FN: " << dumpee.N_FN << "\n";
    output << "N_C:  " << dumpee.N_C  << "\n";
    output << "N_S:  " << dumpee.N_S  << "\n\n";

    output << "b: " << fixnum(dumpee.b) << "\n";
    output << "E: " << fixnum(dumpee.E) << "\n";
    output << "f: " << fixnum(dumpee.f) << "\n\n";

    output << "x:            " << dumpee.x            << "\n";
    output << "z:            " << dumpee.z            << "\n";
    output << "rho:          " << dumpee.rho          << "\n";
    output << "rho_star:     " << dumpee.rho_star     << "\n";
    output << "v:            " << dumpee.v            << "\n";
    output << "h:            " << dumpee.h            << "\n";
    output << "gamma:        " << dumpee.gamma        << "\n";
    output << "gamma_star:   " << dumpee.gamma_star   << "\n";
    output << "mu:           " << dumpee.mu           << "\n";
    output << "mu_star:      " << dumpee.mu_star      << "\n";
    output << "contype:      " << dumpee.contype      << "\n";
    output << "ident:        " << dumpee.ident        << "\n";
    output << "free_ones:    " << dumpee.free_ones    << "\n";
    output << "alpha:        " << dumpee.alpha        << "\n";
    output << "tau:          " << dumpee.tau          << "\n";
    output << "m:            " << dumpee.m            << "\n";
    output << "minv:         " << dumpee.minv         << "\n";
    output << "e_tau:        " << dumpee.e_tau        << "\n";
    output << "e_free_pivot: " << dumpee.e_free_pivot << "\n\n";

    output << "beth_tau: " << dumpee.beth_tau << "\n\n";

    output << "Cache size (MB):     " << dumpee.cache_memsize    << "\n";
    output << "Min cache dimension: " << dumpee.cache_min_rowdim << "\n\n";
    output << "d2c buffer size:     " << dumpee.d2c_buff_size    << "\n\n";

    output << "Kernel: " << dumpee.K << "\n\n";

    FN_EXIT_POINT output;
}

std::ostream &operator<<(std::ostream &output, SVdata_short_summary &dumpee)
{
    FN_ENTRY_POINT;

    output << "SVM data type: S\n\n";

    output << "Flags: " << dumpee.svflags << "\n\n";

    output << "N:    " << dumpee.N    << "\n";
    output << "N_Z:  " << dumpee.N_Z  << "\n";
    output << "N_L:  " << dumpee.N_L  << "\n";
    output << "N_U:  " << dumpee.N_U  << "\n";
    output << "N_F:  " << dumpee.N_F  << "\n";
    output << "N_FP: " << dumpee.N_FP << "\n";
    output << "N_FN: " << dumpee.N_FN << "\n";
    output << "N_C:  " << dumpee.N_C  << "\n";
    output << "N_S:  " << dumpee.N_S  << "\n\n";

    output << "b: " << fixnum(dumpee.b) << "\n";
    output << "E: " << fixnum(dumpee.E) << "\n";
    output << "f: " << fixnum(dumpee.f) << "\n\n";

    output << "contype:      " << dumpee.contype      << "\n";
    output << "alpha:        " << dumpee.alpha        << "\n";
    output << "tau:          " << dumpee.tau          << "\n";
    output << "m:            " << dumpee.m            << "\n";
    output << "minv:         " << dumpee.minv         << "\n";
    output << "e_tau:        " << dumpee.e_tau        << "\n";
    output << "e_free_pivot: " << dumpee.e_free_pivot << "\n\n";

    output << "beth_tau: " << dumpee.beth_tau << "\n\n";

    output << "Cache size (MB):     " << dumpee.cache_memsize    << "\n";
    output << "Min cache dimension: " << dumpee.cache_min_rowdim << "\n\n";
    output << "d2c buffer size:     " << dumpee.d2c_buff_size    << "\n\n";

    FN_EXIT_POINT output;
}



std::istream &operator>>(std::istream &input,  SVdata &source)
{
    FN_ENTRY_POINT;

    char data_type;

    wait_dummy howzat;

    input >> howzat; input >> data_type;

    input >> howzat; input >> source.fix_bias;
    input >> howzat; input >> source.svflags;

    input >> howzat; input >> source.N;
    input >> howzat; input >> source.N_Z;
    input >> howzat; input >> source.N_L;
    input >> howzat; input >> source.N_U;
    input >> howzat; input >> source.N_F;
    input >> howzat; input >> source.N_FP;
    input >> howzat; input >> source.N_FN;
    input >> howzat; input >> source.N_C;
    input >> howzat; input >> source.N_S;

    input >> howzat; input >> source.b;
    input >> howzat; input >> source.E;
    input >> howzat; input >> source.f;

    input >> howzat; input >> source.x;
    input >> howzat; input >> source.z;
    input >> howzat; input >> source.rho;
    input >> howzat; input >> source.rho_star;
    input >> howzat; input >> source.v;
    input >> howzat; input >> source.h;
    input >> howzat; input >> source.gamma;
    input >> howzat; input >> source.gamma_star;
    input >> howzat; input >> source.mu;
    input >> howzat; input >> source.mu_star;
    input >> howzat; input >> source.contype;
    input >> howzat; input >> source.ident;
    input >> howzat; input >> source.free_ones;
    input >> howzat; input >> source.alpha;
    input >> howzat; input >> source.tau;
    input >> howzat; input >> source.m;
    input >> howzat; input >> source.minv;
    input >> howzat; input >> source.e_tau;
    input >> howzat; input >> source.e_free_pivot;

    input >> howzat; input >> source.beth_tau;

    input >> howzat; input >> source.cache_memsize;
    input >> howzat; input >> source.cache_min_rowdim;
    input >> howzat; input >> source.d2c_buff_size;

    input >> howzat; input >> source.K;

    switch ( data_type )
    {
        case 'N':
        {
            input >> howzat; input >> source.G_tau;

            break;
        }

        case 'S':
        {
            /*
               Remake (source.G_tau)
            */

            if ( KERN_CACHE_DYNA(source.svflags) )
            {
                (source.kerncache).setmemsize(&source,source.cache_memsize,source.cache_min_rowdim);
            }

            source.rewrite_G_tau(1);

            break;
        }

        default:
        {
            L_THROW(0);

            break;
        }
    }

    FN_EXIT_POINT input;
}


std::istream &operator>>(std::istream &input,  SVdata_summary &source)
{
    FN_ENTRY_POINT;

    char data_type;

    wait_dummy howzat;

    input >> howzat; input >> data_type;

    input >> howzat; input >> source.fix_bias;
    input >> howzat; input >> source.svflags;

    input >> howzat; input >> source.N;
    input >> howzat; input >> source.N_Z;
    input >> howzat; input >> source.N_L;
    input >> howzat; input >> source.N_U;
    input >> howzat; input >> source.N_F;
    input >> howzat; input >> source.N_FP;
    input >> howzat; input >> source.N_FN;
    input >> howzat; input >> source.N_C;
    input >> howzat; input >> source.N_S;

    input >> howzat; input >> source.b;
    input >> howzat; input >> source.E;
    input >> howzat; input >> source.f;

    input >> howzat; input >> source.x;
    input >> howzat; input >> source.z;
    input >> howzat; input >> source.rho;
    input >> howzat; input >> source.rho_star;
    input >> howzat; input >> source.v;
    input >> howzat; input >> source.h;
    input >> howzat; input >> source.gamma;
    input >> howzat; input >> source.gamma_star;
    input >> howzat; input >> source.mu;
    input >> howzat; input >> source.mu_star;
    input >> howzat; input >> source.contype;
    input >> howzat; input >> source.ident;
    input >> howzat; input >> source.free_ones;
    input >> howzat; input >> source.alpha;
    input >> howzat; input >> source.tau;
    input >> howzat; input >> source.m;
    input >> howzat; input >> source.minv;
    input >> howzat; input >> source.e_tau;
    input >> howzat; input >> source.e_free_pivot;

    input >> howzat; input >> source.beth_tau;

    input >> howzat; input >> source.cache_memsize;
    input >> howzat; input >> source.cache_min_rowdim;
    input >> howzat; input >> source.d2c_buff_size;

    input >> howzat; input >> source.K;

    switch ( data_type )
    {
        case 'N':
        {
            /*
               Clear the input stream, but disgard the result.
            */

            fMATRIX temp;

            input >> howzat; input >> temp;

            break;
        }

        case 'S':
        {
            break;
        }

        default:
        {
            L_THROW(0);

            break;
        }
    }

    FN_EXIT_POINT input;
}




/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/














































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

#define KERNGETROWOPT(__i) (kerncache.opt_getrow(this,((__i)+1)))

void SVdata::enter_opt(void)
{
    FN_ENTRY_POINT;

    long i;

    /*
       Note - this is crude, and makes many assumptions - but on the plus
       side, it's also rather fast.
    */

    fast_z          = &(z[0]);
    fast_rho        = &(rho[0]);
    fast_rho_star   = &(rho_star[0]);
    fast_v          = &(v[0]);
    fast_h          = &(h[0]);
    fast_gamma      = &(gamma[0]);
    fast_gamma_star = &(gamma_star[0]);
    fast_mu         = &(mu[0]);
    fast_mu_star    = &(mu_star[0]);
    fast_alpha      = &(alpha[0]);
    fast_e_tau      = &(e_tau[0]);
    fast_contype    = &(contype[0]);
    fast_tau        = &(tau[0]);
    fast_m          = &(m[0]);
    fast_minv       = &(minv[0]);

    REPNEWB(fast_x,fVECTOR *,N);
    REPNEWB(fast_opt_x,L_DOUBLE *,N);

    for ( i = 0 ; i < N ; i++ )
    {
        fast_x[i]     = &(x[i]);
        fast_opt_x[i] = &((x[i])[0]);
    }

    if ( KERN_CACHE_DYNA(svflags) )
    {
        kerncache.enter_opt(this);
    }

    fast_kernel_fn = K.get_fast_kern(&fast_uc,&fast_uv,(double **) (&fast_covw));

    if ( D2CSMO_OPTIM(svflags) )
    {
        if ( SVM_HINT_PATTERN(svflags) )
        {
            #ifdef CHECK_OPTS
            std::cerr << "d2csmo optimisation\n";
            #endif

            if ( fast_kernel_fn == NULL )
            {
                opt_sel__get_e                   = &SVdata::opt_d2csmo_pattern_get_e;
                opt_sel__get_e_pos               = &SVdata::opt_d2csmo_pattern_get_e_pos;
                opt_sel__get_e_neg               = &SVdata::opt_d2csmo_pattern_get_e_neg;
                opt_sel__get_G_tau               = &SVdata::opt_d2csmo_pattern_get_G_tau;
                opt_sel__get_G_extras            = &SVdata::opt_d2csmo_pattern_get_G_extras;
                opt_sel__get_G_tau_pivot         = &SVdata::opt_d2csmo_pattern_get_G_tau_pivot;
                opt_sel__calc_G_tau_from_scratch = &SVdata::opt_d2csmo_pattern_calc_G_tau_from_scratch;
                opt_sel__step_none               = &SVdata::opt_d2csmo_pattern_step_none;
                opt_sel__step_one                = &SVdata::opt_d2csmo_pattern_step_one;
                opt_sel__step_two                = &SVdata::opt_d2csmo_pattern_step_two;
                opt_sel__step_general_pivot      = &SVdata::opt_d2csmo_pattern_step_general_pivot;
                opt_sel__step_general_pivotx     = &SVdata::opt_d2csmo_pattern_step_general_pivotx;
                opt_sel__constrain_Z             = &SVdata::opt_d2csmo_pattern_constrain_Z;
                opt_sel__constrain_L             = &SVdata::opt_d2csmo_pattern_constrain_L;
                opt_sel__constrain_U             = &SVdata::opt_d2csmo_pattern_constrain_U;
                opt_sel__free_L                  = &SVdata::opt_d2csmo_pattern_free_L;
                opt_sel__free_U                  = &SVdata::opt_d2csmo_pattern_free_U;
                opt_sel__constrain_Z_pivot       = &SVdata::opt_d2csmo_pattern_constrain_Z_pivot;
                opt_sel__constrain_L_pivot       = &SVdata::opt_d2csmo_pattern_constrain_L_pivot;
                opt_sel__constrain_U_pivot       = &SVdata::opt_d2csmo_pattern_constrain_U_pivot;
                opt_sel__free_L_pivot            = &SVdata::opt_d2csmo_pattern_free_L_pivot;
                opt_sel__free_U_pivot            = &SVdata::opt_d2csmo_pattern_free_U_pivot;
            }

            else
            {
                opt_sel__get_e                   = &SVdata::opt_d2csmo_pattern_kfast_get_e;
                opt_sel__get_e_pos               = &SVdata::opt_d2csmo_pattern_kfast_get_e_pos;
                opt_sel__get_e_neg               = &SVdata::opt_d2csmo_pattern_kfast_get_e_neg;
                opt_sel__get_G_tau               = &SVdata::opt_d2csmo_pattern_kfast_get_G_tau;
                opt_sel__get_G_extras            = &SVdata::opt_d2csmo_pattern_kfast_get_G_extras;
                opt_sel__get_G_tau_pivot         = &SVdata::opt_d2csmo_pattern_kfast_get_G_tau_pivot;
                opt_sel__calc_G_tau_from_scratch = &SVdata::opt_d2csmo_pattern_kfast_calc_G_tau_from_scratch;
                opt_sel__step_none               = &SVdata::opt_d2csmo_pattern_kfast_step_none;
                opt_sel__step_one                = &SVdata::opt_d2csmo_pattern_kfast_step_one;
                opt_sel__step_two                = &SVdata::opt_d2csmo_pattern_kfast_step_two;
                opt_sel__step_general_pivot      = &SVdata::opt_d2csmo_pattern_kfast_step_general_pivot;
                opt_sel__step_general_pivotx     = &SVdata::opt_d2csmo_pattern_kfast_step_general_pivotx;
                opt_sel__constrain_Z             = &SVdata::opt_d2csmo_pattern_kfast_constrain_Z;
                opt_sel__constrain_L             = &SVdata::opt_d2csmo_pattern_kfast_constrain_L;
                opt_sel__constrain_U             = &SVdata::opt_d2csmo_pattern_kfast_constrain_U;
                opt_sel__free_L                  = &SVdata::opt_d2csmo_pattern_kfast_free_L;
                opt_sel__free_U                  = &SVdata::opt_d2csmo_pattern_kfast_free_U;
                opt_sel__constrain_Z_pivot       = &SVdata::opt_d2csmo_pattern_kfast_constrain_Z_pivot;
                opt_sel__constrain_L_pivot       = &SVdata::opt_d2csmo_pattern_kfast_constrain_L_pivot;
                opt_sel__constrain_U_pivot       = &SVdata::opt_d2csmo_pattern_kfast_constrain_U_pivot;
                opt_sel__free_L_pivot            = &SVdata::opt_d2csmo_pattern_kfast_free_L_pivot;
                opt_sel__free_U_pivot            = &SVdata::opt_d2csmo_pattern_kfast_free_U_pivot;
            }
        }

        else
        {
            #ifdef CHECK_OPTS
            std::cerr << "d2csmo optimisation\n";
            #endif

            opt_sel__get_e                   = &SVdata::opt_d2csmo_get_e;
            opt_sel__get_e_pos               = &SVdata::opt_d2csmo_get_e_pos;
            opt_sel__get_e_neg               = &SVdata::opt_d2csmo_get_e_neg;
            opt_sel__get_G_tau               = &SVdata::opt_d2csmo_get_G_tau;
            opt_sel__get_G_extras            = &SVdata::opt_d2csmo_get_G_extras;
            opt_sel__get_G_tau_pivot         = &SVdata::opt_d2csmo_get_G_tau_pivot;
            opt_sel__calc_G_tau_from_scratch = &SVdata::opt_d2csmo_calc_G_tau_from_scratch;
            opt_sel__step_none               = &SVdata::opt_d2csmo_step_none;
            opt_sel__step_one                = &SVdata::opt_d2csmo_step_one;
            opt_sel__step_two                = &SVdata::opt_d2csmo_step_two;
            opt_sel__step_general_pivot      = &SVdata::opt_d2csmo_step_general_pivot;
            opt_sel__step_general_pivotx     = &SVdata::opt_d2csmo_step_general_pivotx;
            opt_sel__constrain_Z             = &SVdata::opt_d2csmo_constrain_Z;
            opt_sel__constrain_L             = &SVdata::opt_d2csmo_constrain_L;
            opt_sel__constrain_U             = &SVdata::opt_d2csmo_constrain_U;
            opt_sel__free_L                  = &SVdata::opt_d2csmo_free_L;
            opt_sel__free_U                  = &SVdata::opt_d2csmo_free_U;
            opt_sel__constrain_Z_pivot       = &SVdata::opt_d2csmo_constrain_Z_pivot;
            opt_sel__constrain_L_pivot       = &SVdata::opt_d2csmo_constrain_L_pivot;
            opt_sel__constrain_U_pivot       = &SVdata::opt_d2csmo_constrain_U_pivot;
            opt_sel__free_L_pivot            = &SVdata::opt_d2csmo_free_L_pivot;
            opt_sel__free_U_pivot            = &SVdata::opt_d2csmo_free_U_pivot;
        }
    }

    else if ( ACTIVE_SMALL_OPTIM(svflags) )
    {
        #ifdef CHECK_OPTS
        std::cerr << "small active set optimisation\n";
        #endif

        opt_sel__get_e                   = &SVdata::opt_actsmall_get_e;
        opt_sel__get_e_pos               = &SVdata::opt_actsmall_get_e_pos;
        opt_sel__get_e_neg               = &SVdata::opt_actsmall_get_e_neg;
        opt_sel__get_G_tau               = &SVdata::opt_actsmall_get_G_tau;
        opt_sel__get_G_extras            = &SVdata::opt_actsmall_get_G_extras;
        opt_sel__get_G_tau_pivot         = &SVdata::opt_actsmall_get_G_tau_pivot;
        opt_sel__calc_G_tau_from_scratch = &SVdata::opt_actsmall_calc_G_tau_from_scratch;
        opt_sel__step_none               = &SVdata::opt_actsmall_step_none;
        opt_sel__step_one                = &SVdata::opt_actsmall_step_one;
        opt_sel__step_two                = &SVdata::opt_actsmall_step_two;
        opt_sel__step_general_pivot      = &SVdata::opt_actsmall_step_general_pivot;
        opt_sel__step_general_pivotx     = &SVdata::opt_actsmall_step_general_pivotx;
        opt_sel__constrain_Z             = &SVdata::opt_actsmall_constrain_Z;
        opt_sel__constrain_L             = &SVdata::opt_actsmall_constrain_L;
        opt_sel__constrain_U             = &SVdata::opt_actsmall_constrain_U;
        opt_sel__free_L                  = &SVdata::opt_actsmall_free_L;
        opt_sel__free_U                  = &SVdata::opt_actsmall_free_U;
        opt_sel__constrain_Z_pivot       = &SVdata::opt_actsmall_constrain_Z_pivot;
        opt_sel__constrain_L_pivot       = &SVdata::opt_actsmall_constrain_L_pivot;
        opt_sel__constrain_U_pivot       = &SVdata::opt_actsmall_constrain_U_pivot;
        opt_sel__free_L_pivot            = &SVdata::opt_actsmall_free_L_pivot;
        opt_sel__free_U_pivot            = &SVdata::opt_actsmall_free_U_pivot;
    }

    else if ( ACTIVE_MEDIUM_OPTIM(svflags) )
    {
        #ifdef CHECK_OPTS
        std::cerr << "medium active set optimisation\n";
        #endif

        opt_sel__get_e                   = &SVdata::opt_actmedium_get_e;
        opt_sel__get_e_pos               = &SVdata::opt_actmedium_get_e_pos;
        opt_sel__get_e_neg               = &SVdata::opt_actmedium_get_e_neg;
        opt_sel__get_G_tau               = &SVdata::opt_actmedium_get_G_tau;
        opt_sel__get_G_extras            = &SVdata::opt_actmedium_get_G_extras;
        opt_sel__get_G_tau_pivot         = &SVdata::opt_actmedium_get_G_tau_pivot;
        opt_sel__calc_G_tau_from_scratch = &SVdata::opt_actmedium_calc_G_tau_from_scratch;
        opt_sel__step_none               = &SVdata::opt_actmedium_step_none;
        opt_sel__step_one                = &SVdata::opt_actmedium_step_one;
        opt_sel__step_two                = &SVdata::opt_actmedium_step_two;
        opt_sel__step_general_pivot      = &SVdata::opt_actmedium_step_general_pivot;
        opt_sel__step_general_pivotx     = &SVdata::opt_actmedium_step_general_pivotx;
        opt_sel__constrain_Z             = &SVdata::opt_actmedium_constrain_Z;
        opt_sel__constrain_L             = &SVdata::opt_actmedium_constrain_L;
        opt_sel__constrain_U             = &SVdata::opt_actmedium_constrain_U;
        opt_sel__free_L                  = &SVdata::opt_actmedium_free_L;
        opt_sel__free_U                  = &SVdata::opt_actmedium_free_U;
        opt_sel__constrain_Z_pivot       = &SVdata::opt_actmedium_constrain_Z_pivot;
        opt_sel__constrain_L_pivot       = &SVdata::opt_actmedium_constrain_L_pivot;
        opt_sel__constrain_U_pivot       = &SVdata::opt_actmedium_constrain_U_pivot;
        opt_sel__free_L_pivot            = &SVdata::opt_actmedium_free_L_pivot;
        opt_sel__free_U_pivot            = &SVdata::opt_actmedium_free_U_pivot;
    }

    else
    {
        #ifdef CHECK_OPTS
        std::cerr << "no optimisation\n";
        #endif

        opt_sel__get_e                   = &SVdata::opt_generic_get_e;
        opt_sel__get_e_pos               = &SVdata::opt_generic_get_e_pos;
        opt_sel__get_e_neg               = &SVdata::opt_generic_get_e_neg;
        opt_sel__get_G_tau               = &SVdata::opt_generic_get_G_tau;
        opt_sel__get_G_extras            = &SVdata::opt_generic_get_G_extras;
        opt_sel__get_G_tau_pivot         = &SVdata::opt_generic_get_G_tau_pivot;
        opt_sel__calc_G_tau_from_scratch = &SVdata::opt_generic_calc_G_tau_from_scratch;
        opt_sel__step_none               = &SVdata::opt_generic_step_none;
        opt_sel__step_one                = &SVdata::opt_generic_step_one;
        opt_sel__step_two                = &SVdata::opt_generic_step_two;
        opt_sel__step_general_pivot      = &SVdata::opt_generic_step_general_pivot;
        opt_sel__step_general_pivotx     = &SVdata::opt_generic_step_general_pivotx;
        opt_sel__constrain_Z             = &SVdata::opt_generic_constrain_Z;
        opt_sel__constrain_L             = &SVdata::opt_generic_constrain_L;
        opt_sel__constrain_U             = &SVdata::opt_generic_constrain_U;
        opt_sel__free_L                  = &SVdata::opt_generic_free_L;
        opt_sel__free_U                  = &SVdata::opt_generic_free_U;
        opt_sel__constrain_Z_pivot       = &SVdata::opt_generic_constrain_Z_pivot;
        opt_sel__constrain_L_pivot       = &SVdata::opt_generic_constrain_L_pivot;
        opt_sel__constrain_U_pivot       = &SVdata::opt_generic_constrain_U_pivot;
        opt_sel__free_L_pivot            = &SVdata::opt_generic_free_L_pivot;
        opt_sel__free_U_pivot            = &SVdata::opt_generic_free_U_pivot;
    }

    FN_EXIT_POINT;
}

void SVdata::exit_opt(void)
{
    FN_ENTRY_POINT;

    long i,i_pivot;

    REPDELB(fast_x);
    REPDELB(fast_opt_x);
    fast_x     = NULL;
    fast_opt_x = NULL;

    if ( KERN_CACHE_DYNA(svflags) )
    {
        kerncache.exit_opt();
    }

    if ( D2CSMO_OPTIM(svflags) )
    {
        /*
           Need to fix all counters
        */

        N_Z  = 0;
        N_L  = 0;
        N_U  = 0;
        N_F  = 0;
        N_FP = 0;
        N_FN = 0;
        N_C  = 0;
        N_S  = 0;

        for ( i = 0 ; i < N ; i++ )
        {
            switch ( tau[i] )
            {
                case 0:  { N_Z++;  N_C++;        break; }
                case -1: { N_FN++; N_F++; N_S++; break; }
                case +1: { N_FP++; N_F++; N_S++; break; }
                case -2: { N_L++;  N_C++; N_S++; break; }
                case +2: { N_U++;  N_C++; N_S++; break; }

                default:
                {
                    L_THROW(0);
                    break;
                }
            }
        }
    }

    if ( !SVM_HINT_EFREE_USED(svflags) )
    {
        /*
           Need to make sure that e_free_pivot is the right size at least
        */

        while ( e_free_pivot.get_real_size() < N_F )
        {
            e_free_pivot.addend(0.0);
        }

        while ( e_free_pivot.get_real_size() > N_F )
        {
            e_free_pivot.remove(1);
        }
    }

    if ( !SVM_HINT_USE_PIVOTS(svflags) )
    {
        /*
           Need to fix order vector m, and its inverse m_inv
        */

        i_pivot = 0;

        if ( N_Z > 0 )
        {
            for ( i = 0 ; i < N ; i++ )
            {
                if ( fast_tau[i] == 0 )
                {
                    fast_m[i_pivot] = i+1;

                    i_pivot++;
                }
            }
        }

        if ( N_L > 0 )
        {
            for ( i = 0 ; i < N ; i++ )
            {
                if ( fast_tau[i] == -2 )
                {
                    fast_m[i_pivot] = i+1;

                    i_pivot++;
                }
            }
        }

        if ( N_U > 0 )
        {
            for ( i = 0 ; i < N ; i++ )
            {
                if ( fast_tau[i] == +2 )
                {
                    fast_m[i_pivot] = i+1;

                    i_pivot++;
                }
            }
        }

        if ( N_F > 0 )
        {
            for ( i = 0 ; i < N ; i++ )
            {
                if ( ( fast_tau[i] == -1 ) || ( fast_tau[i] == +1 ) )
                {
                    fast_m[i_pivot] = i+1;

                    i_pivot++;
                }
            }
        }

        for ( i_pivot = 0 ; i_pivot < N ; i_pivot++ )
        {
            fast_minv[fast_m[i_pivot]-1] = i_pivot+1;
        }
    }

    FN_EXIT_POINT;
}

L_DOUBLE SVdata::opt_get_e(long iz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__get_e)(iz);
}

L_DOUBLE SVdata::opt_get_e_pos(long iz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__get_e_pos)(iz);
}

L_DOUBLE SVdata::opt_get_e_neg(long iz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__get_e_neg)(iz);
}

L_DOUBLE SVdata::opt_get_G_tau(long iz, long jz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__get_G_tau)(iz,jz);
}

L_DOUBLE SVdata::opt_get_G_extras(long iz, long jz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__get_G_extras)(iz,jz);
}

L_DOUBLE SVdata::opt_get_G_tau_pivot(long i_pivotz, long j_pivotz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__get_G_tau_pivot)(i_pivotz,j_pivotz);
}

L_DOUBLE SVdata::opt_calc_G_tau_from_scratch(long iz, long jz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__calc_G_tau_from_scratch)(iz,jz);
}

void     SVdata::opt_step_none(L_DOUBLE d_b)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__step_none)(d_b);
}

void     SVdata::opt_step_one(L_DOUBLE d_alpha, long iz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__step_one)(d_alpha,iz);
}

void     SVdata::opt_step_two(L_DOUBLE d_b, L_DOUBLE d_f, L_DOUBLE d_alpha1, L_DOUBLE d_alpha2, long i1z, long i2z)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__step_two)(d_b,d_f,d_alpha1,d_alpha2,i1z,i2z);
}

void     SVdata::opt_step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__step_general_pivot)(start_pivot,finish_pivot,d_alpha_pivot,d_b,d_e_pivot,d_f,n_ignore);
}

void     SVdata::opt_step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__step_general_pivotx)(start_pivot,finish_pivot,d_alpha_pivot,d_b);
}

void     SVdata::opt_constrain_Z(long iz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__constrain_Z)(iz);
}

void     SVdata::opt_constrain_L(long iz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__constrain_L)(iz);
}

void     SVdata::opt_constrain_U(long iz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__constrain_U)(iz);
}

void     SVdata::opt_free_L(long iz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__free_L)(iz);
}

void     SVdata::opt_free_U(long iz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__free_U)(iz);
}

void     SVdata::opt_constrain_Z_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__constrain_Z_pivot)(i_pivotz);
}

void     SVdata::opt_constrain_L_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__constrain_L_pivot)(i_pivotz);
}

void     SVdata::opt_constrain_U_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__constrain_U_pivot)(i_pivotz);
}

void     SVdata::opt_free_L_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__free_L_pivot)(i_pivotz);
}

void     SVdata::opt_free_U_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;
    FN_EXIT_POINT (*this.*opt_sel__free_U_pivot)(i_pivotz);
}



fVECTOR &SVdata::opt_get_e_free_pivot(void)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( SVM_HINT_EFREE_USED(svflags) );

    if ( !E_CACHE_FULL(svflags) && !E_CACHE_SUPP(svflags) && !E_CACHE_FREE(svflags) )
    {
        if ( N_C < N )
        {
            long i_pivotz;

            for ( i_pivotz = N_C ; i_pivotz < N ; i_pivotz++ )
            {
                e_free_pivot[i_pivotz-N_C] = opt_get_e_pivot(i_pivotz);
            }
        }
    }

    FN_EXIT_POINT e_free_pivot;
}

L_DOUBLE SVdata::opt_get_e_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= 0 );
    THROW_ASSERT( i_pivotz <  N );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );

    L_DOUBLE result;

    result = opt_get_e(fast_m[i_pivotz]-1);

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_get_e_neg_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= 0 );
    THROW_ASSERT( i_pivotz <  N );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] <= 0 );

    L_DOUBLE result;

    result = opt_get_e_neg(fast_m[i_pivotz]-1);

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_get_e_pos_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= 0 );
    THROW_ASSERT( i_pivotz <  N );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] >= 0 );

    L_DOUBLE result;

    result = opt_get_e_pos(fast_m[i_pivotz]-1);

    FN_EXIT_POINT result;
}


L_DOUBLE SVdata::opt_get_alpha_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= 0 );
    THROW_ASSERT( i_pivotz <  N );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );

    FN_EXIT_POINT fast_alpha[fast_m[i_pivotz]-1];
}

L_DOUBLE SVdata::opt_get_h_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= 0 );
    THROW_ASSERT( i_pivotz <  N );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );

    FN_EXIT_POINT fast_h[fast_m[i_pivotz]-1];
}

L_DOUBLE SVdata::opt_get_v_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= 0 );
    THROW_ASSERT( i_pivotz <  N );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );

    FN_EXIT_POINT fast_v[fast_m[i_pivotz]-1];
}

int SVdata::opt_get_tau_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= 0 );
    THROW_ASSERT( i_pivotz <  N );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );

    FN_EXIT_POINT fast_tau[fast_m[i_pivotz]-1];
}

int SVdata::opt_get_contype_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= 0 );
    THROW_ASSERT( i_pivotz <  N );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );

    FN_EXIT_POINT fast_contype[fast_m[i_pivotz]-1];
}

long SVdata::opt_get_nbad_factor(void)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) ) );

    long result;

    result = beth_tau.get_nbad();

    FN_EXIT_POINT result;
}

void SVdata::opt_minverse_factor(fVECTOR &d_alpha_pivot, L_DOUBLE &d_b, long z_start, long z_end, int f_zero)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) ) );

    MINVERSE_BETH(d_alpha_pivot,d_b,e_free_pivot,f,z_start,z_end,f_zero);

    FN_EXIT_POINT;
}

void SVdata::opt_minverse_factor(fVECTOR &d_alpha_pivot, long z_start, long z_end)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) ) );

    MINVERSE_BETH_SHORT(d_alpha_pivot,e_free_pivot,z_start,z_end);

    FN_EXIT_POINT;
}

void SVdata::opt_near_invert_factor(fVECTOR &d_alpha_pivot, L_DOUBLE &d_b)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) ) );

    NEARINVERSE_BETH(d_alpha_pivot,d_b);

    FN_EXIT_POINT;
}

void SVdata::opt_near_invert_factor(fVECTOR &d_alpha_pivot)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) ) );

    NEARINVERSE_BETH_SHORT(d_alpha_pivot);

    FN_EXIT_POINT;
}

void SVdata::opt_fix_e_free_pivot(void)
{
    FN_ENTRY_POINT;

    long j_pivotz;

    if ( N_C < N )
    {
        for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
        {
            e_free_pivot[j_pivotz-N_C] = opt_get_e(fast_m[j_pivotz]-1);
        }
    }

    FN_EXIT_POINT;
}



/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

#ifdef ROWPLUS
#undef ROWPLUS
#endif
#define ROWPLUS(_isrow_,_iz_,_jz_) ((_isrow_)[(_iz_)]+opt_generic_get_G_extras(_iz_,_jz_))

L_DOUBLE SVdata::opt_generic_get_e(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );

    L_DOUBLE result;
    long i_pivotz,jz,j_pivotz;

    if ( E_CACHE_FULL(svflags) )
    {
        result = fast_e_tau[iz];

        if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] == 2 ) )
        {
            result -= fast_rho[iz];
        }

        else if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] == 1 ) )
        {
            result += fast_rho_star[iz];
        }
    }

    else if ( E_CACHE_SUPP(svflags) )
    {
        if ( fast_tau[iz] == 0 )
        {
            goto calc_e_from_scratch;
        }

        result = fast_e_tau[iz];
    }

    else if ( E_CACHE_FREE(svflags) )
    {
        if ( !( ( fast_tau[iz] == +1 ) || ( fast_tau[iz] == -1 ) ) )
        {
            goto calc_e_from_scratch;
        }

        result = fast_e_tau[iz];
    }

    else
    {
        calc_e_from_scratch:

        result = b - fast_z[iz];

        if ( N_Z < N )
        {
            if ( SVM_HINT_USE_PIVOTS(svflags) )
            {
                if ( KERN_CACHE_FULL(svflags) )
                {
                    for ( j_pivotz = N_Z ; j_pivotz < N ; j_pivotz++ )
                    {
                        result += G_tau(iz,fast_m[j_pivotz]-1) * fast_alpha[fast_m[j_pivotz]-1];
                    }
                }

                else if ( KERN_CACHE_FREE(svflags) )
                {
                    i_pivotz = fast_minv[iz]-1;

                    if ( i_pivotz < N_C )
                    {
                        for ( j_pivotz = N_Z ; j_pivotz < N ; j_pivotz++ )
                        {
                            result += opt_generic_get_G_tau(iz,fast_m[j_pivotz]-1) * fast_alpha[fast_m[j_pivotz]-1];
                        }
                    }

                    else
                    {
                        if ( N_Z < N_C )
                        {
                            for ( j_pivotz = N_Z ; j_pivotz < N_C ; j_pivotz++ )
                            {
                                result += opt_generic_get_G_tau(iz,fast_m[j_pivotz]-1) * fast_alpha[fast_m[j_pivotz]-1];
                            }
                        }

                        if ( N_C < N )
                        {
                            for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
                            {
                                result += G_tau(i_pivotz-N_C,j_pivotz-N_C) * fast_alpha[fast_m[j_pivotz]-1];
                            }
                        }
                    }
                }

                else
                {
                    for ( j_pivotz = N_Z ; j_pivotz < N ; j_pivotz++ )
                    {
                        result += opt_generic_get_G_tau(iz,fast_m[j_pivotz]-1) * fast_alpha[fast_m[j_pivotz]-1];
                    }
                }
            }

            else
            {
                /*
                   All (non dynamic) kernel caching methods use pivots, so
                   clearly there is no static cache here.  However, dynamic
                   caching methods have optimised alternative functions, so
                   KISS.
                */

                for ( jz = 0 ; jz < N ; jz++ )
                {
                    if ( fast_tau[jz] )
                    {
                        result += opt_generic_get_G_tau(iz,jz) * fast_alpha[jz];
                    }
                }
            }
        }

        if ( !SVM_HINT_RHO_ZERO(svflags) )
        {
            if ( fast_tau[iz] > 0 )
            {
                result += fast_rho[iz];
            }

            else if ( fast_tau[iz] < 0 )
            {
                result -= fast_rho_star[iz];
            }
        }
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_generic_get_e_pos(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] >= 0 );

    L_DOUBLE result;

    result = opt_generic_get_e(iz);

    if ( fast_tau[iz] == 0 )
    {
        if ( !SVM_HINT_RHO_ZERO(svflags) )
        {
            result += fast_rho[iz];
        }

        if ( !SVM_HINT_MU_ZERO(svflags) )
        {
            result += E*fast_mu[iz];
        }
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_generic_get_e_neg(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] <= 0 );

    L_DOUBLE result;

    result = opt_generic_get_e(iz);

    if ( fast_tau[iz] == 0 )
    {
        if ( !SVM_HINT_RHO_ZERO(svflags) )
        {
            result -= fast_rho_star[iz];
        }

        if ( !SVM_HINT_MU_ZERO(svflags) )
        {
            result -= E*fast_mu_star[iz];
        }
    }

    FN_EXIT_POINT result;
}


L_DOUBLE SVdata::opt_generic_get_G_tau(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    L_DOUBLE result;
    long i_pivotz,j_pivotz;

    if ( KERN_CACHE_FULL(svflags) )
    {
        result = G_tau(iz,jz);
    }

    else if ( KERN_CACHE_FREE(svflags) )
    {
        THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );

        i_pivotz = fast_minv[iz]-1;
        j_pivotz = fast_minv[jz]-1;

        if ( ( i_pivotz >= N_C ) && ( j_pivotz >= N_C ) )
        {
            result = G_tau(i_pivotz-N_C,j_pivotz-N_C);
        }

        else
        {
            result = opt_generic_calc_G_tau_from_scratch(iz,jz);
        }
    }

    else if ( KERN_CACHE_DYNA(svflags) )
    {
        result = (kerncache).opt_getval_statcache(this,iz+1,jz+1);

        result += opt_generic_get_G_extras(iz,jz);
    }

    else
    {
        result = opt_generic_calc_G_tau_from_scratch(iz,jz);
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_generic_get_G_extras(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    L_DOUBLE result;

    result = 0.0;

    if ( !SVM_HINT_GAMMA_ZERO(svflags) )
    {
        if ( iz == jz )
        {
            if ( fast_tau[iz] > 0 )
            {
                result += fast_gamma[iz];
            }

            else if ( fast_tau[iz] < 0 )
            {
                result += fast_gamma_star[iz];
            }
        }
    }

    if ( !SVM_HINT_MU_ZERO(svflags) )
    {
        if ( ( fast_tau[iz] > 0 ) && ( fast_tau[jz] > 0 ) )
        {
            result += fast_mu[iz] * fast_mu[jz];
        }

        else if ( ( fast_tau[iz] > 0 ) && ( fast_tau[jz] < 0 ) )
        {
            result -= fast_mu[iz] * fast_mu_star[jz];
        }

        else if ( ( fast_tau[iz] < 0 ) && ( fast_tau[jz] > 0 ) )
        {
            result -= fast_mu_star[iz] * fast_mu[jz];
        }

        else if ( ( fast_tau[iz] < 0 ) && ( fast_tau[jz] < 0 ) )
        {
            result += fast_mu_star[iz] * fast_mu_star[jz];
        }
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_generic_get_G_tau_pivot(long i_pivotz, long j_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= 0 );
    THROW_ASSERT( i_pivotz <  N );
    THROW_ASSERT( j_pivotz >= 0 );
    THROW_ASSERT( j_pivotz <  N );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );

    L_DOUBLE result;

    if ( KERN_CACHE_FULL(svflags) )
    {
        result = G_tau(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
    }

    else if ( KERN_CACHE_FREE(svflags) )
    {
        if ( ( i_pivotz >= N_C ) && ( j_pivotz >= N_C ) )
        {
            result = G_tau(i_pivotz-N_C,j_pivotz-N_C);
        }

        else
        {
            result = opt_generic_calc_G_tau_from_scratch(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
        }
    }

    else
    {
        result = opt_generic_calc_G_tau_from_scratch(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_generic_calc_G_tau_from_scratch(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    L_DOUBLE result;

    result = K.kernel(*(fast_x[iz]),*(fast_x[jz]),fast_opt_x[iz],fast_opt_x[jz],0.0,0.0,iz+1,jz+1,fast_z[iz],fast_z[jz]) + opt_generic_get_G_extras(iz,jz);

    FN_EXIT_POINT result;
}


void SVdata::opt_generic_step_none(L_DOUBLE d_b)
{
    FN_ENTRY_POINT;

    long jz,j_pivotz,j_pivot_startz;

    if ( SVM_HINT_USE_PIVOTS(svflags) )
    {
        if ( E_CACHE_FULL(svflags) )
        {
            for ( j_pivotz = 0 ; j_pivotz < N ; j_pivotz++ )
            {
                fast_e_tau[fast_m[j_pivotz]-1] += d_b;
            }

            if ( SVM_HINT_EFREE_USED(svflags) )
            {
                opt_fix_e_free_pivot();
            }

            j_pivot_startz = 0;
        }

        else if ( E_CACHE_SUPP(svflags) )
        {
            for ( j_pivotz = N_Z ; j_pivotz < N ; j_pivotz++ )
            {
                fast_e_tau[fast_m[j_pivotz]-1] += d_b;
            }

            if ( SVM_HINT_EFREE_USED(svflags) )
            {
                opt_fix_e_free_pivot();
            }

            j_pivot_startz = N_Z;
        }

        else if ( E_CACHE_FREE(svflags) )
        {
            for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
            {
                fast_e_tau[fast_m[j_pivotz]-1] += d_b;
            }

            if ( SVM_HINT_EFREE_USED(svflags) )
            {
                opt_fix_e_free_pivot();
            }

            j_pivot_startz = N_C;
        }
    }

    else
    {
        if ( E_CACHE_FULL(svflags) )
        {
            for ( jz = 0 ; jz < N ; jz++ )
            {
                fast_e_tau[jz] += d_b;
            }
        }

        else if ( E_CACHE_SUPP(svflags) )
        {
            if ( N_Z < N )
            {
                for ( jz = 0 ; jz < N ; jz++ )
                {
                    if ( fast_tau[jz] )
                    {
                        fast_e_tau[jz] += d_b;
                    }
                }
            }
        }

        else if ( E_CACHE_FREE(svflags) )
        {
            if ( N_C < N )
            {
                for ( jz = 0 ; jz < N ; jz++ )
                {
                    if ( ( fast_tau[jz] == +1 ) || ( fast_tau[jz] == -1 ) )
                    {
                        fast_e_tau[jz] += d_b;
                    }
                }
            }
        }
    }

    b += d_b;

    FN_EXIT_POINT;
}

void SVdata::opt_generic_step_one(L_DOUBLE d_alpha, long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( ( fast_tau[iz] == -1 ) || ( fast_tau[iz] == +1 ) );

    long i_pivotz,jz,j_pivotz,j_pivot_startz;
    L_DOUBLE *isrow;

    fast_alpha[iz] += d_alpha;

    if ( !SVM_HINT_MU_ZERO(svflags) )
    {
        if ( fast_tau[iz] > 0 )
        {
            E += d_alpha * fast_mu[iz];
        }

        else
        {
            E -= d_alpha * fast_mu_star[iz];
        }
    }

    f += d_alpha;

    if ( SVM_HINT_USE_PIVOTS(svflags) )
    {
        if ( E_CACHE_FULL(svflags) )
        {
            j_pivot_startz = 0;
        }

        else if ( E_CACHE_SUPP(svflags) )
        {
            j_pivot_startz = N_Z;
        }

        else if ( E_CACHE_FREE(svflags) )
        {
            j_pivot_startz = N_C;
        }

        else
        {
            j_pivot_startz = N;
        }

        if ( j_pivot_startz < N )
        {
            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j_pivotz = j_pivot_startz ; j_pivotz < N ; j_pivotz++ )
                {
                    fast_e_tau[fast_m[j_pivotz]-1] += G_tau(fast_m[j_pivotz]-1,iz)*d_alpha;
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                /*
                   Note j_pivot_start = 1,N_Z+1,N_C+1
                */

                i_pivotz = fast_minv[iz]-1;

                /*
                   Take care - internal functions can call this for variables
                   not technically "free"
                */

                if ( j_pivot_startz < N )
                {
                    for ( j_pivotz = j_pivot_startz ; j_pivotz < N ; j_pivotz++ )
                    {
                        fast_e_tau[fast_m[j_pivotz]-1] += opt_generic_get_G_tau(fast_m[j_pivotz]-1,iz)*d_alpha;
                    }
                }

                else
                {
                    if ( j_pivot_startz < N_C )
                    {
                        for ( j_pivotz = j_pivot_startz ; j_pivotz < N_C ; j_pivotz++ )
                        {
                            fast_e_tau[fast_m[j_pivotz]-1] += opt_generic_get_G_tau(fast_m[j_pivotz]-1,iz)*d_alpha;
                        }
                    }

                    for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
                    {
                        fast_e_tau[fast_m[j_pivotz]-1] += G_tau(j_pivotz-N_C,i_pivotz-N_C)*d_alpha;
                    }
                }
            }

            else if ( KERN_CACHE_DYNA(svflags) )
            {
                isrow = KERNGETROWOPT(iz);

                for ( j_pivotz = j_pivot_startz ; j_pivotz < N ; j_pivotz++ )
                {
                    fast_e_tau[fast_m[j_pivotz]-1] += ROWPLUS(isrow,fast_m[j_pivotz]-1,iz) * d_alpha;
                }
            }

            else
            {
                for ( j_pivotz = j_pivot_startz ; j_pivotz < N ; j_pivotz++ )
                {
                    fast_e_tau[fast_m[j_pivotz]-1] += opt_generic_get_G_tau(fast_m[j_pivotz]-1,iz)*d_alpha;
                }
            }

            if ( SVM_HINT_EFREE_USED(svflags) )
            {
                opt_fix_e_free_pivot();
            }
        }
    }

    else
    {
        /*
           Must have N_F > 0 here (or otherwise what's being stepped?)
        */

        if ( E_CACHE_FULL(svflags) )
        {
            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( jz = 0 ; jz < N ; jz++ )
                {
                    fast_e_tau[jz] += G_tau(jz,iz) * d_alpha;
                }
            }

            else if ( KERN_CACHE_DYNA(svflags) )
            {
                isrow = KERNGETROWOPT(iz);

                for ( jz = 0 ; jz < N ; jz++ )
                {
                    fast_e_tau[jz] += ROWPLUS(isrow,jz,iz) * d_alpha;
                }
            }

            else
            {
                for ( jz = 0 ; jz < N ; jz++ )
                {
                    fast_e_tau[jz] += opt_generic_get_G_tau(jz,iz) * d_alpha;
                }
            }
        }

        else if ( E_CACHE_SUPP(svflags) )
        {
            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( jz = 0 ; jz < N ; jz++ )
                {
                    if ( fast_tau[jz] )
                    {
                        fast_e_tau[jz] += G_tau(jz,iz) * d_alpha;
                    }
                }
            }

            else if ( KERN_CACHE_DYNA(svflags) )
            {
                isrow = KERNGETROWOPT(iz);

                for ( jz = 0 ; jz < N ; jz++ )
                {
                    if ( fast_tau[jz] )
                    {
                        fast_e_tau[jz] += ROWPLUS(isrow,jz,iz) * d_alpha;
                    }
                }
            }

            else
            {
                for ( jz = 0 ; jz < N ; jz++ )
                {
                    if ( fast_tau[jz] )
                    {
                        fast_e_tau[jz] += opt_generic_get_G_tau(jz,iz) * d_alpha;
                    }
                }
            }
        }

        else if ( E_CACHE_FREE(svflags) )
        {
            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( jz = 0 ; jz < N ; jz++ )
                {
                    if ( ( fast_tau[jz] == +1 ) || ( fast_tau[jz] == -1 ) )
                    {
                        fast_e_tau[jz] += G_tau(jz,iz) * d_alpha;
                    }
                }
            }

            else if ( KERN_CACHE_DYNA(svflags) )
            {
                isrow = KERNGETROWOPT(iz);

                for ( jz = 0 ; jz < N ; jz++ )
                {
                    if ( ( fast_tau[jz] == +1 ) || ( fast_tau[jz] == -1 ) )
                    {
                        fast_e_tau[jz] += ROWPLUS(isrow,jz,iz) * d_alpha;
                    }
                }
            }

            else
            {
                for ( jz = 0 ; jz < N ; jz++ )
                {
                    if ( ( fast_tau[jz] == +1 ) || ( fast_tau[jz] == -1 ) )
                    {
                        fast_e_tau[jz] += opt_generic_get_G_tau(jz,iz) * d_alpha;
                    }
                }
            }
        }
    }

    FN_EXIT_POINT;
}

void SVdata::opt_generic_step_two(L_DOUBLE d_b, L_DOUBLE d_f, L_DOUBLE d_alpha1, L_DOUBLE d_alpha2, long i1z, long i2z)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i1z >= 0 );
    THROW_ASSERT( i1z <  N );
    THROW_ASSERT( i2z >= 0 );
    THROW_ASSERT( i2z <  N );
    THROW_ASSERT( ( fast_tau[i1z] == -1 ) || ( fast_tau[i1z] == +1 ) );
    THROW_ASSERT( ( fast_tau[i2z] == -1 ) || ( fast_tau[i2z] == +1 ) );

    long i1_pivotz;
    long i2_pivotz;
    long jz,j_pivotz,j_pivot_startz;
    L_DOUBLE *isrow1;
    L_DOUBLE *isrow2;

    fast_alpha[i1z] += d_alpha1;
    fast_alpha[i2z] += d_alpha2;

    f += d_f;

    if ( !SVM_HINT_MU_ZERO(svflags) )
    {
        if ( fast_tau[i1z] > 0 )
        {
            E += d_alpha1 * fast_mu[i1z];
        }

        else
        {
            E -= d_alpha1 * fast_mu_star[i1z];
        }

        if ( fast_tau[i2z] > 0 )
        {
            E += d_alpha2 * fast_mu[i2z];
        }

        else
        {
            E -= d_alpha2 * fast_mu_star[i2z];
        }
    }

    if ( SVM_HINT_USE_PIVOTS(svflags) )
    {
        if ( E_CACHE_FULL(svflags) )
        {
            j_pivot_startz = 0;
        }

        else if ( E_CACHE_SUPP(svflags) )
        {
            j_pivot_startz = N_Z;
        }

        else if ( E_CACHE_FREE(svflags) )
        {
            j_pivot_startz = N_C;
        }

        else
        {
            j_pivot_startz = N;
        }

        if ( j_pivot_startz < N )
        {
            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j_pivotz = j_pivot_startz ; j_pivotz < N ; j_pivotz++ )
                {
                    fast_e_tau[fast_m[j_pivotz]-1] += G_tau(fast_m[j_pivotz]-1,i1z)*d_alpha1;
                    fast_e_tau[fast_m[j_pivotz]-1] += G_tau(fast_m[j_pivotz]-1,i2z)*d_alpha2;
                    fast_e_tau[fast_m[j_pivotz]-1] += d_b;
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                i1_pivotz = fast_minv[i1z]-1;
                i2_pivotz = fast_minv[i2z]-1;

                /*
                   Note j_pivot_start = 1,N_Z+1,N_C+1
                */

                if ( j_pivot_startz < N_C )
                {
                    for ( j_pivotz = j_pivot_startz ; j_pivotz < N_C ; j_pivotz++ )
                    {
                        fast_e_tau[fast_m[j_pivotz]-1] += opt_generic_get_G_tau(fast_m[j_pivotz]-1,i1z)*d_alpha1;
                        fast_e_tau[fast_m[j_pivotz]-1] += opt_generic_get_G_tau(fast_m[j_pivotz]-1,i2z)*d_alpha2;
                        fast_e_tau[fast_m[j_pivotz]-1] += d_b;
                    }
                }

                for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
                {
                    fast_e_tau[fast_m[j_pivotz]-1] += G_tau(j_pivotz-N_C,i1_pivotz-N_C)*d_alpha1;
                    fast_e_tau[fast_m[j_pivotz]-1] += G_tau(j_pivotz-N_C,i2_pivotz-N_C)*d_alpha2;
                    fast_e_tau[fast_m[j_pivotz]-1] += d_b;
                }
            }

            else if ( KERN_CACHE_DYNA(svflags) )
            {
                isrow1 = KERNGETROWOPT(i1z);
                isrow2 = KERNGETROWOPT(i2z);

                for ( j_pivotz = j_pivot_startz ; j_pivotz < N ; j_pivotz++ )
                {
                    fast_e_tau[fast_m[j_pivotz]-1] += ROWPLUS(isrow1,fast_m[j_pivotz]-1,i1z)*d_alpha1;
                    fast_e_tau[fast_m[j_pivotz]-1] += ROWPLUS(isrow2,fast_m[j_pivotz]-1,i2z)*d_alpha2;
                    fast_e_tau[fast_m[j_pivotz]-1] += d_b;
                }
            }

            else
            {
                for ( j_pivotz = j_pivot_startz ; j_pivotz < N ; j_pivotz++ )
                {
                    fast_e_tau[fast_m[j_pivotz]-1] += opt_generic_get_G_tau(fast_m[j_pivotz]-1,i1z)*d_alpha1;
                    fast_e_tau[fast_m[j_pivotz]-1] += opt_generic_get_G_tau(fast_m[j_pivotz]-1,i2z)*d_alpha2;
                    fast_e_tau[fast_m[j_pivotz]-1] += d_b;
                }
            }

            if ( SVM_HINT_EFREE_USED(svflags) )
            {
                opt_fix_e_free_pivot();
            }
        }
    }

    else
    {
        if ( E_CACHE_FULL(svflags) )
        {
            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( jz = 0 ; jz < N ; jz++ )
                {
                    fast_e_tau[jz] += G_tau(jz,i1z) * d_alpha1;
                    fast_e_tau[jz] += G_tau(jz,i2z) * d_alpha2;
                    fast_e_tau[jz] += d_b;
                }
            }

            else if ( KERN_CACHE_DYNA(svflags) )
            {
                isrow1 = KERNGETROWOPT(i1z);
                isrow2 = KERNGETROWOPT(i2z);

                for ( jz = 0 ; jz < N ; jz++ )
                {
                    fast_e_tau[jz] += ROWPLUS(isrow1,jz,i1z) * d_alpha1;
                    fast_e_tau[jz] += ROWPLUS(isrow2,jz,i2z) * d_alpha2;
                    fast_e_tau[jz] += d_b;
                }
            }

            else
            {
                for ( jz = 0 ; jz < N ; jz++ )
                {
                    fast_e_tau[jz] += opt_generic_get_G_tau(jz,i1z) * d_alpha1;
                    fast_e_tau[jz] += opt_generic_get_G_tau(jz,i2z) * d_alpha2;
                    fast_e_tau[jz] += d_b;
                }
            }
        }

        else if ( E_CACHE_SUPP(svflags) )
        {
            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( jz = 0 ; jz < N ; jz++ )
                {
                    if ( fast_tau[jz] )
                    {
                        fast_e_tau[jz] += G_tau(jz,i1z) * d_alpha1;
                        fast_e_tau[jz] += G_tau(jz,i2z) * d_alpha2;
                        fast_e_tau[jz] += d_b;
                    }
                }
            }

            else if ( KERN_CACHE_DYNA(svflags) )
            {
                isrow1 = KERNGETROWOPT(i1z);
                isrow2 = KERNGETROWOPT(i2z);

                for ( jz = 0 ; jz < N ; jz++ )
                {
                    if ( fast_tau[jz] )
                    {
                        fast_e_tau[jz] += ROWPLUS(isrow1,jz,i1z) * d_alpha1;
                        fast_e_tau[jz] += ROWPLUS(isrow2,jz,i2z) * d_alpha2;
                        fast_e_tau[jz] += d_b;
                    }
                }
            }

            else
            {
                for ( jz = 0 ; jz < N ; jz++ )
                {
                    if ( fast_tau[jz] )
                    {
                        fast_e_tau[jz] += opt_generic_get_G_tau(jz,i1z) * d_alpha1;
                        fast_e_tau[jz] += opt_generic_get_G_tau(jz,i2z) * d_alpha2;
                        fast_e_tau[jz] += d_b;
                    }
                }
            }
        }

        else if ( E_CACHE_FREE(svflags) )
        {
            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( jz = 0 ; jz < N ; jz++ )
                {
                    if ( ( fast_tau[jz] == +1 ) || ( fast_tau[jz] == -1 ) )
                    {
                        fast_e_tau[jz] += G_tau(jz,i1z) * d_alpha1;
                        fast_e_tau[jz] += G_tau(jz,i2z) * d_alpha2;
                        fast_e_tau[jz] += d_b;
                    }
                }
            }

            else if ( KERN_CACHE_DYNA(svflags) )
            {
                isrow1 = KERNGETROWOPT(i1z);
                isrow2 = KERNGETROWOPT(i2z);

                for ( jz = 0 ; jz < N ; jz++ )
                {
                    if ( ( fast_tau[jz] == +1 ) || ( fast_tau[jz] == -1 ) )
                    {
                        fast_e_tau[jz] += ROWPLUS(isrow1,jz,i1z) * d_alpha1;
                        fast_e_tau[jz] += ROWPLUS(isrow2,jz,i2z) * d_alpha2;
                        fast_e_tau[jz] += d_b;
                    }
                }
            }

            else
            {
                for ( jz = 0 ; jz < N ; jz++ )
                {
                    if ( ( fast_tau[jz] == +1 ) || ( fast_tau[jz] == -1 ) )
                    {
                        fast_e_tau[jz] += opt_generic_get_G_tau(jz,i1z) * d_alpha1;
                        fast_e_tau[jz] += opt_generic_get_G_tau(jz,i2z) * d_alpha2;
                        fast_e_tau[jz] += d_b;
                    }
                }
            }
        }
    }

    b += d_b;

    FN_EXIT_POINT;
}

void SVdata::opt_generic_step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore)
{
    /*
       but see also other instance of this function.
    */

    FN_ENTRY_POINT;

    THROW_ASSERT( d_alpha_pivot.get_effective_size() == d_e_pivot.get_effective_size() );
    THROW_ASSERT( d_alpha_pivot.get_effective_size() == finish_pivot-start_pivot+1 );
    THROW_ASSERT( start_pivot >= N_C+1 );
    THROW_ASSERT( start_pivot <= N );
    THROW_ASSERT( finish_pivot >= N_C+1 );
    THROW_ASSERT( finish_pivot <= N );
    THROW_ASSERT( start_pivot <= finish_pivot );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );

    long i_pivot,j_pivot,i_pivot_start;

    for ( i_pivot = start_pivot ; i_pivot <= finish_pivot ; i_pivot++ )
    {
        fast_alpha[fast_m[i_pivot-1]-1] += d_alpha_pivot[i_pivot-start_pivot];
    }

    if ( E_CACHE_FULL(svflags) || E_CACHE_SUPP(svflags) || E_CACHE_FREE(svflags) )
    {
        for ( i_pivot = start_pivot+n_ignore ; i_pivot <= finish_pivot ; i_pivot++ )
        {
            fast_e_tau[fast_m[i_pivot-1]-1] += d_e_pivot[i_pivot-start_pivot];
        }
    }

    /*
       update E and f

       E = mu_(star)^T |alpha|
       f = 1^T alpha
    */

    f += d_f;

    if ( !SVM_HINT_MU_ZERO(svflags) )
    {
        for ( i_pivot = start_pivot ; i_pivot <= finish_pivot ; i_pivot++ )
        {
            if ( fast_tau[fast_m[i_pivot-1]-1] > 0 )
            {
                E += d_alpha_pivot[i_pivot-start_pivot] * fast_mu[fast_m[i_pivot-1]-1];
            }

            else
            {
                E -= d_alpha_pivot[i_pivot-start_pivot] * fast_mu_star[fast_m[i_pivot-1]-1];
            }
        }
    }

    /*
       update e_tau
    */

    if ( E_CACHE_FULL(svflags) )
    {
        i_pivot_start = 1;
    }

    else if ( E_CACHE_SUPP(svflags) )
    {
        i_pivot_start = N_Z+1;
    }

    else if ( E_CACHE_FREE(svflags) )
    {
        i_pivot_start = N_C+1;
    }

    else
    {
        i_pivot_start = N+1;
    }

    if ( ( start_pivot > 1 ) && ( i_pivot_start < start_pivot ) )
    {
        for ( i_pivot = i_pivot_start ; i_pivot < start_pivot ; i_pivot++ )
        {
            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    fast_e_tau[fast_m[i_pivot-1]-1] += G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    fast_e_tau[fast_m[i_pivot-1]-1] += opt_generic_get_G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else if ( KERN_CACHE_DYNA(svflags) )
            {
                /*
                   This cannot happen with the current set of optimisation
                   schemes.
                */

                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    fast_e_tau[fast_m[i_pivot-1]-1] += opt_generic_get_G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else
            {
                /*
                   This cannot happen with the current set of optimisation
                   schemes.
                */

                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    fast_e_tau[fast_m[i_pivot-1]-1] += opt_generic_get_G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            fast_e_tau[fast_m[i_pivot-1]-1] += d_b;
        }
    }

    if ( i_pivot_start < finish_pivot+1 )
    {
        i_pivot_start = finish_pivot+1;
    }

    if ( i_pivot_start <= N )
    {
        for ( i_pivot = i_pivot_start ; i_pivot <= N ; i_pivot++ )
        {
            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    fast_e_tau[fast_m[i_pivot-1]-1] += G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    fast_e_tau[fast_m[i_pivot-1]-1] += opt_generic_get_G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else if ( KERN_CACHE_DYNA(svflags) )
            {
                /*
                   This cannot happen with the current set of optimisation
                   schemes.
                */

                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    fast_e_tau[fast_m[i_pivot-1]-1] += opt_generic_get_G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else
            {
                /*
                   This cannot happen with the current set of optimisation
                   schemes.
                */

                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    fast_e_tau[fast_m[i_pivot-1]-1] += opt_generic_get_G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            fast_e_tau[fast_m[i_pivot-1]-1] += d_b;
        }
    }

    b += d_b;

    if ( SVM_HINT_EFREE_USED(svflags) )
    {
        opt_fix_e_free_pivot();
    }

    FN_EXIT_POINT;
}

void SVdata::opt_generic_step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b)
{
    /*
       but see also other instance of this function.
    */

    FN_ENTRY_POINT;

    THROW_ASSERT( d_alpha_pivot.get_effective_size() == finish_pivot-start_pivot+1 );
    THROW_ASSERT( start_pivot >= N_C+1 );
    THROW_ASSERT( start_pivot <= N );
    THROW_ASSERT( finish_pivot >= N_C+1 );
    THROW_ASSERT( finish_pivot <= N );
    THROW_ASSERT( start_pivot <= finish_pivot );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );

    long i_pivot,j_pivot,i_pivot_start;

    for ( i_pivot = start_pivot ; i_pivot <= finish_pivot ; i_pivot++ )
    {
        fast_alpha[fast_m[i_pivot-1]-1] += d_alpha_pivot[i_pivot-start_pivot];
    }

    /*
       update E and f

       E = mu_(star)^T |alpha|
       f = 1^T alpha
    */

    if ( !SVM_HINT_MU_ZERO(svflags) )
    {
        for ( i_pivot = start_pivot ; i_pivot <= finish_pivot ; i_pivot++ )
        {
            if ( fast_tau[fast_m[i_pivot-1]-1] > 0 )
            {
                E += d_alpha_pivot[i_pivot-start_pivot] * fast_mu[fast_m[i_pivot-1]-1];
            }

            else
            {
                E -= d_alpha_pivot[i_pivot-start_pivot] * fast_mu_star[fast_m[i_pivot-1]-1];
            }
        }
    }

    for ( i_pivot = start_pivot ; i_pivot <= finish_pivot ; i_pivot++ )
    {
        f += d_alpha_pivot[i_pivot-start_pivot];
    }

    /*
       update e_tau
    */

    if ( E_CACHE_FULL(svflags) )
    {
        i_pivot_start = 1;
    }

    else if ( E_CACHE_SUPP(svflags) )
    {
        i_pivot_start = N_Z+1;
    }

    else if ( E_CACHE_FREE(svflags) )
    {
        i_pivot_start = N_C+1;
    }

    else
    {
        i_pivot_start = N+1;
    }

    if ( i_pivot_start <= N )
    {
        for ( i_pivot = i_pivot_start ; i_pivot <= N ; i_pivot++ )
        {
            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    fast_e_tau[fast_m[i_pivot-1]-1] += G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    fast_e_tau[fast_m[i_pivot-1]-1] += opt_generic_get_G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else if ( KERN_CACHE_DYNA(svflags) )
            {
                /*
                   This cannot happen with the current set of optimisation
                   schemes.
                */

                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    fast_e_tau[fast_m[i_pivot-1]-1] += opt_generic_get_G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            else
            {
                /*
                   This cannot happen with the current set of optimisation
                   schemes.
                */

                for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
                {
                    fast_e_tau[fast_m[i_pivot-1]-1] += opt_generic_get_G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
                }
            }

            fast_e_tau[fast_m[i_pivot-1]-1] += d_b;
        }

        if ( SVM_HINT_EFREE_USED(svflags) )
        {
            opt_fix_e_free_pivot();
        }
    }

    b += d_b;

    FN_EXIT_POINT;
}



/*
   Pivotting.
*/

void SVdata::opt_generic_constrain_Z(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != -2 );
    THROW_ASSERT( fast_tau[iz] != +2 );

    long jz;

    if ( SVM_HINT_USE_PIVOTS(svflags) )
    {
        opt_generic_constrain_Z_pivot(fast_minv[iz]-1);
    }

    else if ( fast_tau[iz] != 0 )
    {
        if ( fast_tau[iz] == -1 )
        {
            /*
               Fixup G_tau - this operation takes advantage of the symmetry
               property of the hessian matrix class, so that the operation
               must be done for only one of the row/column.
            */

            if ( KERN_CACHE_FULL(svflags) )
            {
                if ( !SVM_HINT_GAMMA_ZERO(svflags) )
                {
                    G_tau(iz,iz) -= fast_gamma_star[iz];
                }

                if ( !SVM_HINT_MU_ZERO(svflags) )
                {
                    for ( jz = 0 ; jz < N ; jz++ )
                    {
                        if ( fast_tau[jz] < 0 )
                        {
                            G_tau(iz,jz) -= ( fast_mu_star[iz] * fast_mu_star[jz] );
                        }

                        else if ( fast_tau[jz] > 0 )
                        {
                            G_tau(iz,jz) += ( fast_mu_star[iz] * fast_mu[jz] );
                        }
                    }
                }
            }

            /*
               Fixup e_tau
            */

            if ( E_CACHE_FULL(svflags) && !SVM_HINT_RHO_ZERO(svflags) && ( fast_contype[iz] != 1 ) )
            {
                fast_e_tau[iz] += fast_rho_star[iz];
            }

            if ( E_CACHE_FULL(svflags) && !SVM_HINT_MU_ZERO(svflags) )
            {
                fast_e_tau[iz] += E * fast_mu_star[iz];
            }

            N_FN--;
        }

        else
        {
            /*
               Fixup G_tau - this operation takes advantage of the symmetry
               property of the hessian matrix class, so that the operation
               must be done for only one of the row/column.
            */

            if ( KERN_CACHE_FULL(svflags) )
            {
                if ( !SVM_HINT_GAMMA_ZERO(svflags) )
                {
                    G_tau(iz,iz) -= fast_gamma[iz];
                }

                if ( !SVM_HINT_MU_ZERO(svflags) )
                {
                    for ( jz = 0 ; jz < N ; jz++ )
                    {
                        if ( fast_tau[jz] < 0 )
                        {
                            G_tau(iz,jz) += ( fast_mu[iz] * fast_mu_star[jz] );
                        }

                        else if ( fast_tau[jz] > 0 )
                        {
                            G_tau(iz,jz) -= ( fast_mu[iz] * fast_mu[jz] );
                        }
                    }
                }
            }

            /*
               Fixup e_tau
            */

            if ( E_CACHE_FULL(svflags) && !SVM_HINT_RHO_ZERO(svflags) && ( fast_contype[iz] != 2 ) )
            {
                fast_e_tau[iz] -= fast_rho[iz];
            }

            if ( E_CACHE_FULL(svflags) && !SVM_HINT_MU_ZERO(svflags) )
            {
                fast_e_tau[iz] -= E * fast_mu[iz];
            }

            N_FP--;
        }

        /*
           Fixup tau
        */

        fast_tau[iz] = 0;

        N_Z++;
        N_F--;
        N_C++;
        N_S--;
    }

    FN_EXIT_POINT;
}

void SVdata::opt_generic_constrain_L(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != 0 );
    THROW_ASSERT( fast_tau[iz] != +1 );
    THROW_ASSERT( fast_tau[iz] != +2 );

    if ( SVM_HINT_USE_PIVOTS(svflags) )
    {
        opt_generic_constrain_L_pivot(fast_minv[iz]-1);
    }

    else if ( fast_tau[iz] != -2 )
    {
        /*
           Neither G_tau nor e_tau will be changed by this operation.
        */

        /*
           Fixup tau
        */

        fast_tau[iz] = -2;

        N_L++;
        N_F--;
        N_FN--;
        N_C++;
    }

    FN_EXIT_POINT;
}

void SVdata::opt_generic_constrain_U(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != 0 );
    THROW_ASSERT( fast_tau[iz] != -1 );
    THROW_ASSERT( fast_tau[iz] != -2 );

    if ( SVM_HINT_USE_PIVOTS(svflags) )
    {
        opt_generic_constrain_U_pivot(fast_minv[iz]-1);
    }

    else if ( fast_tau[iz] != +2 )
    {
        /*
           Neither G_tau nor e_tau will be changed by this operation.
        */

        /*
           Fixup tau
        */

        fast_tau[iz] = +2;

        N_U++;
        N_F--;
        N_FP--;
        N_C++;
    }

    FN_EXIT_POINT;
}

void SVdata::opt_generic_free_L(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != +1 );
    THROW_ASSERT( fast_tau[iz] != +2 );
    THROW_ASSERT( fast_contype[iz] !=  0 );
    THROW_ASSERT( fast_contype[iz] != +2 );

    long jz;
    int old_tau;

    if ( SVM_HINT_USE_PIVOTS(svflags) )
    {
        opt_generic_free_L_pivot(fast_minv[iz]);
    }

    else if ( fast_tau[iz] != -1 )
    {
        old_tau = fast_tau[iz];

        /*
           If e is entering the cache, will need to ensure that it has a
           valid value before we start.
        */

        if ( ( E_CACHE_SUPP(svflags) && ( old_tau == 0 ) ) || ( E_CACHE_FREE(svflags) ) )
        {
            fast_e_tau[iz] = opt_generic_get_e(iz);

            if ( ( fast_contype[iz] == 1 ) && !SVM_HINT_RHO_ZERO(svflags) )
            {
                fast_e_tau[iz] -= fast_rho_star[iz];
            }
        }

        if ( old_tau == 0 )
        {
            N_Z--;
            N_F++;
            N_FN++;
            N_C--;
            N_S++;
        }

        else
        {
            N_L--;
            N_F++;
            N_FN++;
            N_C--;
        }

        /*
           Fixup tau
        */

        fast_tau[iz] = -1;

        if ( old_tau == 0 )
        {
            /*
               Fixup G_tau - this operation takes advantage of the symmetry
               property of the hessian matrix class, so that the operation
               must be done for only one of the row/column.
            */

            if ( KERN_CACHE_FULL(svflags) )
            {
                if ( !SVM_HINT_GAMMA_ZERO(svflags) )
                {
                    G_tau(iz,iz) += fast_gamma_star[iz];
                }

                if ( !SVM_HINT_MU_ZERO(svflags) )
                {
                    for ( jz = 0 ; jz < N ; jz++ )
                    {
                        if ( fast_tau[jz] < 0 )
                        {
                            G_tau(iz,jz) += ( fast_mu_star[iz] * fast_mu_star[jz] );
                        }

                        else if ( fast_tau[jz] > 0 )
                        {
                            G_tau(iz,jz) -= ( fast_mu_star[iz] * fast_mu[jz] );
                        }
                    }
                }
            }

            /*
               Fixup e_tau
            */

            if ( ( E_CACHE_FULL(svflags) || E_CACHE_SUPP(svflags) || E_CACHE_FREE(svflags) ) && !SVM_HINT_RHO_ZERO(svflags) && ( fast_contype[iz] != 1 ) )
            {
                fast_e_tau[iz] -= fast_rho_star[iz];
            }

            if ( ( E_CACHE_FULL(svflags) || E_CACHE_SUPP(svflags) || E_CACHE_FREE(svflags) ) && !SVM_HINT_MU_ZERO(svflags) )
            {
                fast_e_tau[iz] -= E * fast_mu_star[iz];
            }
        }
    }

    FN_EXIT_POINT;
}

void SVdata::opt_generic_free_U(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != -1 );
    THROW_ASSERT( fast_tau[iz] != -2 );
    THROW_ASSERT( fast_contype[iz] !=  0 );
    THROW_ASSERT( fast_contype[iz] != +1 );

    long jz;
    int old_tau;

    if ( SVM_HINT_USE_PIVOTS(svflags) )
    {
        opt_generic_free_U_pivot(fast_minv[iz]);
    }

    else if ( fast_tau[iz] != +1 )
    {
        old_tau = fast_tau[iz];

        /*
           If e is entering the cache, will need to ensure that it has a
           valid value before we start.
        */

        if ( ( E_CACHE_SUPP(svflags) && ( old_tau == 0 ) ) || ( E_CACHE_FREE(svflags) ) )
        {
            fast_e_tau[iz] = opt_generic_get_e(iz);

            if ( ( fast_contype[iz] == 2 ) && !SVM_HINT_RHO_ZERO(svflags) )
            {
                fast_e_tau[iz] += fast_rho[iz];
            }
        }

        if ( old_tau == 0 )
        {
            N_Z--;
            N_F++;
            N_FP++;
            N_C--;
            N_S++;
        }

        else
        {
            N_U--;
            N_F++;
            N_FP++;
            N_C--;
        }

        /*
           Fixup tau
        */

        fast_tau[iz] = +1;

        if ( old_tau == 0 )
        {
            /*
               Fixup G_tau - this operation takes advantage of the symmetry
               property of the hessian matrix class, so that the operation
               must be done for only one of the row/column.
            */

            if ( KERN_CACHE_FULL(svflags) )
            {
                if ( !SVM_HINT_GAMMA_ZERO(svflags) )
                {
                    G_tau(iz,iz) += fast_gamma[iz];
                }

                if ( !SVM_HINT_MU_ZERO(svflags) )
                {
                    for ( jz = 0 ; jz < N ; jz++ )
                    {
                        if ( fast_tau[jz] < 0 )
                        {
                            G_tau(iz,jz) -= ( fast_mu[iz] * fast_mu_star[jz] );
                        }

                        else if ( fast_tau[jz] > 0 )
                        {
                            G_tau(iz,jz) += ( fast_mu[iz] * fast_mu[jz] );
                        }
                    }
                }
            }

            /*
               Fixup e_tau
            */

            if ( ( E_CACHE_FULL(svflags) || E_CACHE_SUPP(svflags) || E_CACHE_FREE(svflags) ) && !SVM_HINT_RHO_ZERO(svflags) && ( fast_contype[iz] != 2 ) )
            {
                fast_e_tau[iz] += fast_rho[iz];
            }

            if ( ( E_CACHE_FULL(svflags) || E_CACHE_SUPP(svflags) || E_CACHE_FREE(svflags) ) && !SVM_HINT_MU_ZERO(svflags) )
            {
                fast_e_tau[iz] += E * fast_mu[iz];
            }
        }
    }

    FN_EXIT_POINT;
}

void SVdata::opt_generic_constrain_Z_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= N_C );
    THROW_ASSERT( i_pivotz <  N   );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != -2 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != +2 );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );

    long iz,jz,j_pivotz;

    if ( fast_tau[fast_m[i_pivotz]-1] != 0 )
    {
        if ( SVM_HINT_EFREE_USED(svflags) )
        {
            /*
               fix e_free_pivot
            */

            e_free_pivot.remove(i_pivotz+1-N_C);
        }

        iz = fast_m[i_pivotz]-1;

        if ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) )
        {
            /*
               Fixup beth_tau
            */

            SHRINK_BETH(i_pivotz+1-N_C);
        }

        if ( KERN_CACHE_FREE(svflags) )
        {
            G_tau.remove(i_pivotz+1-N_C);
        }

        if ( fast_tau[iz] == -1 )
        {
            /*
               Fixup G_tau - this operation takes advantage of the symmetry
               property of the hessian matrix class, so that the operation
               must be done for only one of the row/column.
            */

            if ( KERN_CACHE_FULL(svflags) )
            {
                if ( !SVM_HINT_GAMMA_ZERO(svflags) )
                {
                    G_tau(iz,iz) -= fast_gamma_star[iz];
                }

                if ( !SVM_HINT_MU_ZERO(svflags) )
                {
                    for ( j_pivotz = N_Z ; j_pivotz < N ; j_pivotz++ )
                    {
                        jz = fast_m[j_pivotz]-1;

                        if ( fast_tau[jz] < 0 )
                        {
                            G_tau(iz,jz) -= ( fast_mu_star[iz] * fast_mu_star[jz] );
                        }

                        else
                        {
                            G_tau(iz,jz) += ( fast_mu_star[iz] * fast_mu[jz] );
                        }
                    }
                }
            }

            /*
               Fixup e_tau
            */

            if ( E_CACHE_FULL(svflags) && !SVM_HINT_RHO_ZERO(svflags) && ( fast_contype[iz] != 1 ) )
            {
                fast_e_tau[iz] += fast_rho_star[iz];
            }

            if ( E_CACHE_FULL(svflags) && !SVM_HINT_MU_ZERO(svflags) )
            {
                fast_e_tau[iz] += E * fast_mu_star[iz];
            }

            N_FN--;
        }

        else
        {
            /*
               Fixup G_tau - this operation takes advantage of the symmetry
               property of the hessian matrix class, so that the operation
               must be done for only one of the row/column.
            */

            if ( KERN_CACHE_FULL(svflags) )
            {
                if ( !SVM_HINT_GAMMA_ZERO(svflags) )
                {
                    G_tau(iz,iz) -= fast_gamma[iz];
                }

                if ( !SVM_HINT_MU_ZERO(svflags) )
                {
                    for ( j_pivotz = N_Z ; j_pivotz < N ; j_pivotz++ )
                    {
                        jz = fast_m[j_pivotz]-1;

                        if ( fast_tau[jz] < 0 )
                        {
                            G_tau(iz,jz) += ( fast_mu[iz] * fast_mu_star[jz] );
                        }

                        else
                        {
                            G_tau(iz,jz) -= ( fast_mu[iz] * fast_mu[jz] );
                        }
                    }
                }
            }

            /*
               Fixup e_tau
            */

            if ( E_CACHE_FULL(svflags) && !SVM_HINT_RHO_ZERO(svflags) && ( fast_contype[iz] != 2 ) )
            {
                fast_e_tau[iz] -= fast_rho[iz];
            }

            if ( E_CACHE_FULL(svflags) && !SVM_HINT_MU_ZERO(svflags) )
            {
                fast_e_tau[iz] -= E * fast_mu[iz];
            }

            N_FP--;
        }

        /*
           Fixup tau
        */

        fast_tau[iz] = 0;

        /*
           Fixup the order vector
        */

        m.bswap(N_Z+1,i_pivotz+1);

        for ( j_pivotz = N_Z ; j_pivotz <= i_pivotz ; j_pivotz++ )
        {
            fast_minv[fast_m[j_pivotz]-1] = j_pivotz+1;
        }

        N_Z++;
        N_F--;
        N_C++;
        N_S--;
    }

    FN_EXIT_POINT;
}

void SVdata::opt_generic_constrain_L_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= N_C );
    THROW_ASSERT( i_pivotz <  N   );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != 0 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != +1 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != +2 );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );

    long iz,j_pivotz;

    if ( fast_tau[fast_m[i_pivotz]-1] != -2 )
    {
        if ( SVM_HINT_EFREE_USED(svflags) )
        {
            /*
               fix e_free_pivot
            */

            e_free_pivot.remove(i_pivotz+1-N_C);
        }

        iz = fast_m[i_pivotz]-1;

        if ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) )
        {
            /*
               Fixup beth_tau
            */

            SHRINK_BETH(i_pivotz+1-N_C);
        }

        if ( KERN_CACHE_FREE(svflags) )
        {
            G_tau.remove(i_pivotz+1-N_C);
        }

        /*
           Neither G_tau nor e_tau will be changed by this operation.
        */

        /*
           Fixup tau
        */

        fast_tau[iz] = -2;

        /*
           Fixup the order vector
        */

        m.bswap(N_Z+N_L+1,i_pivotz+1);

        for ( j_pivotz = N_Z+N_L ; j_pivotz <= i_pivotz ; j_pivotz++ )
        {
            fast_minv[fast_m[j_pivotz]-1] = j_pivotz+1;
        }

        N_L++;
        N_F--;
        N_FN--;
        N_C++;
    }

    FN_EXIT_POINT;
}

void SVdata::opt_generic_constrain_U_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= N_C );
    THROW_ASSERT( i_pivotz <  N   );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != -2 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != -1 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != 0 );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );

    long iz,j_pivotz;

    if ( fast_tau[fast_m[i_pivotz]-1] != +2 )
    {
        if ( SVM_HINT_EFREE_USED(svflags) )
        {
            /*
               fix e_free_pivot
            */

            e_free_pivot.remove(i_pivotz+1-N_C);
        }

        iz = fast_m[i_pivotz]-1;

        if ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) )
        {
            /*
               Fixup beth_tau
            */

            SHRINK_BETH(i_pivotz+1-N_C);
        }

        if ( KERN_CACHE_FREE(svflags) )
        {
            G_tau.remove(i_pivotz+1-N_C);
        }

        /*
           Neither G_tau nor e_tau will be changed by this operation.
        */

        /*
           Fixup tau
        */

        fast_tau[iz] = +2;

        /*
           Fixup the order vector
        */

        m.bswap(N_C+1,i_pivotz+1);

        for ( j_pivotz = N_C ; j_pivotz <= i_pivotz ; j_pivotz++ )
        {
            fast_minv[fast_m[j_pivotz]-1] = j_pivotz;
        }

        N_U++;
        N_F--;
        N_FP--;
        N_C++;
    }

    FN_EXIT_POINT;
}

void SVdata::opt_generic_free_L_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= 0 );
    THROW_ASSERT( i_pivotz <= N_Z+N_L );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != 1 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != 2 );
    THROW_ASSERT( fast_contype[fast_m[i_pivotz]-1] != 0 );
    THROW_ASSERT( fast_contype[fast_m[i_pivotz]-1] != 2 );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );

    long iz,jz,j_pivotz;
    int old_tau;

    if ( fast_tau[fast_m[i_pivotz]-1] != -1 )
    {
        iz = fast_m[i_pivotz]-1;

        old_tau = fast_tau[iz];

        /*
           If e is entering the cache, will need to ensure that it has a
           valid value before we start.
        */

        if ( ( E_CACHE_SUPP(svflags) && ( old_tau == 0 ) ) || ( E_CACHE_FREE(svflags) ) )
        {
            fast_e_tau[iz] = opt_generic_get_e(iz);

            if ( fast_contype[iz] == 1 )
            {
                fast_e_tau[iz] -= fast_rho_star[iz];
            }
        }

        /*
           Pivotting will be needed first.
        */

        if ( old_tau == 0 )
        {
            N_Z--;
            N_F++;
            N_FN++;
            N_C--;
            N_S++;
        }

        else
        {
            N_L--;
            N_F++;
            N_FN++;
            N_C--;
        }

        /*
           Fixup tau
        */

        fast_tau[iz] = -1;

        /*
           Fixup the order vector
        */

        m.fswap(i_pivotz+1,N);

        for ( j_pivotz = i_pivotz ; j_pivotz < N ; j_pivotz++ )
        {
            fast_minv[fast_m[j_pivotz]-1] = j_pivotz+1;
        }

        /*
           At this point, the point has moved, so reset i_pivot.
        */

        i_pivotz = N-1;

        if ( old_tau == 0 )
        {
            /*
               Fixup G_tau - this operation takes advantage of the symmetry
               property of the hessian matrix class, so that the operation
               must be done for only one of the row/column.
            */

            if ( KERN_CACHE_FULL(svflags) )
            {
                if ( !SVM_HINT_GAMMA_ZERO(svflags) )
                {
                    G_tau(iz,iz) += fast_gamma_star[iz];
                }

                if ( !SVM_HINT_MU_ZERO(svflags) )
                {
                    for ( j_pivotz = N_Z ; j_pivotz < N ; j_pivotz++ )
                    {
                        jz = fast_m[j_pivotz]-1;

                        if ( fast_tau[jz] < 0 )
                        {
                            G_tau(iz,jz) += ( fast_mu_star[iz] * fast_mu_star[jz] );
                        }

                        else
                        {
                            G_tau(iz,jz) -= ( fast_mu_star[iz] * fast_mu[jz] );
                        }
                    }
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                fVECTOR ax('x',N_F);

                for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
                {
                    ax[j_pivotz-N_C] = opt_generic_calc_G_tau_from_scratch(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
                }

                G_tau.addend(ax[N_F-1],ax);
            }

            /*
               Fixup e_tau
            */

            if ( ( E_CACHE_FULL(svflags) || E_CACHE_SUPP(svflags) || E_CACHE_FREE(svflags) ) && !SVM_HINT_RHO_ZERO(svflags) && ( fast_contype[iz] != 1 ) )
            {
                fast_e_tau[iz] -= fast_rho_star[iz];
            }

            if ( ( E_CACHE_FULL(svflags) || E_CACHE_SUPP(svflags) || E_CACHE_FREE(svflags) ) && !SVM_HINT_MU_ZERO(svflags) )
            {
                fast_e_tau[iz] -= E * fast_mu_star[iz];
            }
        }

        else if ( KERN_CACHE_FREE(svflags) )
        {
            fVECTOR ax('x',N_F);

            for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
            {
                ax[j_pivotz-N_C] = opt_generic_calc_G_tau_from_scratch(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
            }

            G_tau.addend(ax[N_F-1],ax);
        }

        if ( SVM_HINT_EFREE_USED(svflags) )
        {
            /*
               Add new element to e_free_pivot
            */

            e_free_pivot.addend(opt_generic_get_e(iz));
        }

        if ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) )
        {
            /*
               Fixup beth_tau
            */

            fVECTOR ax('x',N_F);

            /*
               Implicitly, if a factorisation is kept then G_tau must
               be either FULL or FREE
            */

            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
                {
                    ax[j_pivotz-N_C] = G_tau(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
                }
            }

            else
            {
                for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
                {
                    ax[j_pivotz-N_C] = G_tau(i_pivotz-N_C,j_pivotz-N_C);
                }
            }

            ADDEND_BETH(ax);
        }
    }

    FN_EXIT_POINT;
}

void SVdata::opt_generic_free_U_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( ( ( i_pivotz >= 0 ) && ( i_pivotz < N_Z ) ) || ( ( i_pivotz >= N_Z+N_L ) && ( i_pivotz < N_C ) ) );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != -1 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != -2 );
    THROW_ASSERT( fast_contype[fast_m[i_pivotz]-1] != 0 );
    THROW_ASSERT( fast_contype[fast_m[i_pivotz]-1] != 1 );
    THROW_ASSERT( SVM_HINT_USE_PIVOTS(svflags) );

    long iz,jz,j_pivotz;
    int old_tau;

    if ( fast_tau[fast_m[i_pivotz]-1] != +1 )
    {
        iz = fast_m[i_pivotz]-1;

        old_tau = fast_tau[iz];

        /*
           If e is entering the cache, will need to ensure that it has a
           valid value before we start.
        */

        if ( ( E_CACHE_SUPP(svflags) && ( old_tau == 0 ) ) || ( E_CACHE_FREE(svflags) ) )
        {
            fast_e_tau[iz] = opt_generic_get_e(iz);

            if ( fast_contype[iz] == 2 )
            {
                fast_e_tau[iz] += fast_rho[iz];
            }
        }

        if ( old_tau == 0 )
        {
            N_Z--;
            N_F++;
            N_FP++;
            N_C--;
            N_S++;
        }

        else
        {
            N_U--;
            N_F++;
            N_FP++;
            N_C--;
        }

        /*
           Fixup tau
        */

        fast_tau[iz] = +1;

        /*
           Fixup the order vector
        */

        m.fswap(i_pivotz+1,N);

        for ( j_pivotz = i_pivotz ; j_pivotz < N ; j_pivotz++ )
        {
            fast_minv[fast_m[j_pivotz]-1] = j_pivotz+1;
        }

        /*
           At this point, the point has moved, so reset i_pivot.
        */

        i_pivotz = N-1;

        if ( old_tau == 0 )
        {
            /*
               Fixup G_tau - this operation takes advantage of the symmetry
               property of the hessian matrix class, so that the operation
               must be done for only one of the row/column.
            */

            if ( KERN_CACHE_FULL(svflags) )
            {
                if ( !SVM_HINT_GAMMA_ZERO(svflags) )
                {
                    G_tau(iz,iz) += fast_gamma[iz];
                }

                if ( !SVM_HINT_MU_ZERO(svflags) )
                {
                    for ( j_pivotz = N_Z ; j_pivotz < N ; j_pivotz++ )
                    {
                        jz = fast_m[j_pivotz]-1;

                        if ( fast_tau[jz] < 0 )
                        {
                            G_tau(iz,jz) -= ( fast_mu[iz] * fast_mu_star[jz] );
                        }

                        else
                        {
                            G_tau(iz,jz) += ( fast_mu[iz] * fast_mu[jz] );
                        }
                    }
                }
            }

            else if ( KERN_CACHE_FREE(svflags) )
            {
                fVECTOR ax('x',N_F);

                for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
                {
                    ax[j_pivotz-N_C] = opt_generic_calc_G_tau_from_scratch(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
                }

                G_tau.addend(ax[N_F-1],ax);
            }

            /*
               Fixup e_tau
            */

            if ( ( E_CACHE_FULL(svflags) || E_CACHE_SUPP(svflags) || E_CACHE_FREE(svflags) ) && !SVM_HINT_RHO_ZERO(svflags) && ( fast_contype[iz] != 2 ) )
            {
                fast_e_tau[iz] += fast_rho[iz];
            }

            if ( ( E_CACHE_FULL(svflags) || E_CACHE_SUPP(svflags) || E_CACHE_FREE(svflags) ) && !SVM_HINT_MU_ZERO(svflags) )
            {
                fast_e_tau[iz] += E * fast_mu[iz];
            }
        }

        else if ( KERN_CACHE_FREE(svflags) )
        {
            fVECTOR ax('x',N_F);

            for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
            {
                ax[j_pivotz-N_C] = opt_generic_calc_G_tau_from_scratch(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
            }

            G_tau.addend(ax[N_F-1],ax);
        }

        if ( SVM_HINT_EFREE_USED(svflags) )
        {
            /*
               Add new element to e_free_pivot
            */

            e_free_pivot.addend(opt_generic_get_e(iz));
        }

        if ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) )
        {
            /*
               Fixup beth_tau
            */

            fVECTOR ax('x',N_F);

            if ( KERN_CACHE_FULL(svflags) )
            {
                for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
                {
                    ax[j_pivotz-N_C] = G_tau(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
                }
            }

            else
            {
                for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
                {
                    ax[j_pivotz-N_C] = G_tau(i_pivotz-N_C,j_pivotz-N_C);
                }
            }

            ADDEND_BETH(ax);
        }
    }

    FN_EXIT_POINT;
}


/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/


L_DOUBLE SVdata::opt_actsmall_get_e(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );

    L_DOUBLE result;

    result = fast_e_tau[iz];

    if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] == 2 ) )
    {
        result -= fast_rho[iz];
    }

    else if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] == 1 ) )
    {
        result += fast_rho_star[iz];
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_actsmall_get_e_pos(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] >= 0 );

    L_DOUBLE result;

    result = fast_e_tau[iz];

    if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] != 2 ) )
    {
        result += fast_rho[iz];
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_actsmall_get_e_neg(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] <= 0 );

    L_DOUBLE result;

    result = fast_e_tau[iz];

    if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] != 1 ) )
    {
        result -= fast_rho_star[iz];
    }

    FN_EXIT_POINT result;
}


L_DOUBLE SVdata::opt_actsmall_get_G_tau(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    FN_EXIT_POINT G_tau(iz,jz);
}

L_DOUBLE SVdata::opt_actsmall_get_G_extras(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    FN_EXIT_POINT 0.0;

    iz = 1;
    jz = 1;
}

L_DOUBLE SVdata::opt_actsmall_get_G_tau_pivot(long i_pivotz, long j_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= 0 );
    THROW_ASSERT( i_pivotz <  N );
    THROW_ASSERT( j_pivotz >= 0 );
    THROW_ASSERT( j_pivotz <  N );

    FN_EXIT_POINT G_tau(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
}

L_DOUBLE SVdata::opt_actsmall_calc_G_tau_from_scratch(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    L_DOUBLE result;

    result = K.kernel(*(fast_x[iz]),*(fast_x[jz]),fast_opt_x[iz],fast_opt_x[jz],0.0,0.0,iz+1,jz+1,fast_z[iz],fast_z[jz]);

    FN_EXIT_POINT result;
}


void SVdata::opt_actsmall_step_none(L_DOUBLE d_b)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    d_b = 0.0;

    FN_EXIT_POINT;
}

void SVdata::opt_actsmall_step_one(L_DOUBLE d_alpha, long iz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    d_alpha = 0.0;
    iz      = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_actsmall_step_two(L_DOUBLE d_b, L_DOUBLE d_f, L_DOUBLE d_alpha1, L_DOUBLE d_alpha2, long i1z, long i2z)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    d_b      = 0.0;
    d_f      = 0.0;
    d_alpha1 = 0.0;
    d_alpha2 = 0.0;
    i1z      = 0;
    i2z      = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_actsmall_step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( d_alpha_pivot.get_effective_size() == d_e_pivot.get_effective_size() );
    THROW_ASSERT( d_alpha_pivot.get_effective_size() == finish_pivot-start_pivot+1 );
    THROW_ASSERT( start_pivot >= N_C+1 );
    THROW_ASSERT( start_pivot <= N );
    THROW_ASSERT( finish_pivot >= N_C+1 );
    THROW_ASSERT( finish_pivot <= N );
    THROW_ASSERT( start_pivot <= finish_pivot );

    long i_pivot,j_pivot;

    if ( n_ignore > 0 )
    {
        for ( i_pivot = start_pivot ; i_pivot < start_pivot+n_ignore ; i_pivot++ )
        {
            fast_alpha[fast_m[i_pivot-1]-1] += d_alpha_pivot[i_pivot-start_pivot];
        }
    }

    if ( start_pivot+n_ignore <= finish_pivot )
    {
        for ( i_pivot = start_pivot+n_ignore ; i_pivot <= finish_pivot ; i_pivot++ )
        {
            fast_alpha[fast_m[i_pivot-1]-1] += d_alpha_pivot[i_pivot-start_pivot];
            fast_e_tau[fast_m[i_pivot-1]-1] += d_e_pivot[i_pivot-start_pivot];
        }
    }

    f += d_f;

    if ( start_pivot > 1 )
    {
        for ( i_pivot = 1 ; i_pivot < start_pivot ; i_pivot++ )
        {
            for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
            {
                fast_e_tau[fast_m[i_pivot-1]-1] += G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
            }

            fast_e_tau[fast_m[i_pivot-1]-1] += d_b;
        }
    }

    if ( finish_pivot < N )
    {
        for ( i_pivot = finish_pivot+1 ; i_pivot <= N ; i_pivot++ )
        {
            for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
            {
                fast_e_tau[fast_m[i_pivot-1]-1] += G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
            }

            fast_e_tau[fast_m[i_pivot-1]-1] += d_b;
        }
    }

    b += d_b;

    opt_fix_e_free_pivot();

    FN_EXIT_POINT;
}

void SVdata::opt_actsmall_step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( d_alpha_pivot.get_effective_size() == finish_pivot-start_pivot+1 );
    THROW_ASSERT( start_pivot >= N_C+1 );
    THROW_ASSERT( start_pivot <= N );
    THROW_ASSERT( finish_pivot >= N_C+1 );
    THROW_ASSERT( finish_pivot <= N );
    THROW_ASSERT( start_pivot <= finish_pivot );

    long i_pivot,j_pivot;

    for ( i_pivot = start_pivot ; i_pivot <= finish_pivot ; i_pivot++ )
    {
        fast_alpha[fast_m[i_pivot-1]-1] += d_alpha_pivot[i_pivot-start_pivot];
        f += d_alpha_pivot[i_pivot-start_pivot];
    }

    for ( i_pivot = 1 ; i_pivot <= N ; i_pivot++ )
    {
        for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
        {
            fast_e_tau[fast_m[i_pivot-1]-1] += G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
        }

        fast_e_tau[fast_m[i_pivot-1]-1] += d_b;
    }

    b += d_b;

    opt_fix_e_free_pivot();

    FN_EXIT_POINT;
}

void SVdata::opt_actsmall_constrain_Z(long iz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    iz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_actsmall_constrain_L(long iz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    iz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_actsmall_constrain_U(long iz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    iz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_actsmall_free_L(long iz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    iz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_actsmall_free_U(long iz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    iz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_actsmall_constrain_Z_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= N_C );
    THROW_ASSERT( i_pivotz <  N   );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != -2 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != +2 );

    long iz,j_pivotz;

    if ( fast_tau[fast_m[i_pivotz]-1] != 0 )
    {
        e_free_pivot.remove(i_pivotz+1-N_C);

        iz = fast_m[i_pivotz]-1;

        SHRINK_BETH(i_pivotz+1-N_C);

        if ( fast_tau[iz] == -1 )
        {
            if ( fast_contype[iz] != 1 )
            {
                fast_e_tau[iz] += fast_rho_star[iz];
            }

            N_FN--;
        }

        else
        {
            if ( fast_contype[iz] != 2 )
            {
                fast_e_tau[iz] -= fast_rho[iz];
            }

            N_FP--;
        }

        fast_tau[iz] = 0;

        m.bswap(N_Z+1,i_pivotz+1);

        for ( j_pivotz = N_Z ; j_pivotz <= i_pivotz ; j_pivotz++ )
        {
            fast_minv[fast_m[j_pivotz]-1] = j_pivotz+1;
        }

        N_Z++;
        N_F--;
        N_C++;
        N_S--;
    }

    FN_EXIT_POINT;
}

void SVdata::opt_actsmall_constrain_L_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= N_C );
    THROW_ASSERT( i_pivotz <  N   );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != 0 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != +1 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != +2 );

    long iz,j_pivotz;

    if ( fast_tau[fast_m[i_pivotz]-1] != -2 )
    {
        e_free_pivot.remove(i_pivotz+1-N_C);

        iz = fast_m[i_pivotz]-1;

        SHRINK_BETH(i_pivotz+1-N_C);

        fast_tau[iz] = -2;

        m.bswap(N_Z+N_L+1,i_pivotz+1);

        for ( j_pivotz = N_Z+N_L ; j_pivotz <= i_pivotz ; j_pivotz++ )
        {
            fast_minv[fast_m[j_pivotz]-1] = j_pivotz+1;
        }

        N_L++;
        N_F--;
        N_FN--;
        N_C++;
    }

    FN_EXIT_POINT;
}

void SVdata::opt_actsmall_constrain_U_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= N_C );
    THROW_ASSERT( i_pivotz <  N   );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != -2 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != -1 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != 0 );

    long iz,j_pivotz;

    if ( fast_tau[fast_m[i_pivotz]-1] != +2 )
    {
        e_free_pivot.remove(i_pivotz+1-N_C);

        iz = fast_m[i_pivotz]-1;

        SHRINK_BETH(i_pivotz+1-N_C);

        fast_tau[iz] = +2;

        m.bswap(N_C+1,i_pivotz+1);

        for ( j_pivotz = N_C ; j_pivotz <= i_pivotz ; j_pivotz++ )
        {
            fast_minv[fast_m[j_pivotz]-1] = j_pivotz+1;
        }

        N_U++;
        N_F--;
        N_FP--;
        N_C++;
    }

    FN_EXIT_POINT;
}

void SVdata::opt_actsmall_free_L_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= 0 );
    THROW_ASSERT( i_pivotz <= N_Z+N_L );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != 1 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != 2 );
    THROW_ASSERT( fast_contype[fast_m[i_pivotz]-1] != 0 );
    THROW_ASSERT( fast_contype[fast_m[i_pivotz]-1] != 2 );

    long iz,j_pivotz;
    int old_tau;

    if ( fast_tau[fast_m[i_pivotz]-1] != -1 )
    {
        iz = fast_m[i_pivotz]-1;

        old_tau = fast_tau[iz];

        if ( old_tau == 0 )
        {
            N_Z--;
            N_F++;
            N_FN++;
            N_C--;
            N_S++;
        }

        else
        {
            N_L--;
            N_F++;
            N_FN++;
            N_C--;
        }

        fast_tau[iz] = -1;

        m.fswap(i_pivotz+1,N);

        for ( j_pivotz = i_pivotz ; j_pivotz < N ; j_pivotz++ )
        {
            fast_minv[fast_m[j_pivotz]-1] = j_pivotz+1;
        }

        i_pivotz = N-1;

        if ( ( old_tau == 0 ) && ( fast_contype[iz] != 1 ) )
        {
            fast_e_tau[iz] -= fast_rho_star[iz];
        }

        e_free_pivot.addend(fast_e_tau[iz]);

        fVECTOR ax('x',N_F);

        for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
        {
            ax[j_pivotz-N_C] = G_tau(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
        }

        ADDEND_BETH(ax);
    }

    FN_EXIT_POINT;
}

void SVdata::opt_actsmall_free_U_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( ( ( i_pivotz >= 0 ) && ( i_pivotz < N_Z ) ) || ( ( i_pivotz >= N_Z+N_L ) && ( i_pivotz < N_C ) ) );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != -1 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != -2 );
    THROW_ASSERT( fast_contype[fast_m[i_pivotz]-1] != 0 );
    THROW_ASSERT( fast_contype[fast_m[i_pivotz]-1] != 1 );

    long iz,j_pivotz;
    int old_tau;

    if ( fast_tau[fast_m[i_pivotz]-1] != +1 )
    {
        iz = fast_m[i_pivotz]-1;

        old_tau = fast_tau[iz];

        if ( old_tau == 0 )
        {
            N_Z--;
            N_F++;
            N_FP++;
            N_C--;
            N_S++;
        }

        else
        {
            N_U--;
            N_F++;
            N_FP++;
            N_C--;
        }

        fast_tau[iz] = +1;

        m.fswap(i_pivotz+1,N);

        for ( j_pivotz = i_pivotz ; j_pivotz < N ; j_pivotz++ )
        {
            fast_minv[fast_m[j_pivotz]-1] = j_pivotz+1;
        }

        i_pivotz = N-1;

        if ( ( old_tau == 0 ) && ( fast_contype[iz] != 2 ) )
        {
            fast_e_tau[iz] += fast_rho[iz];
        }

        e_free_pivot.addend(fast_e_tau[iz]);

        fVECTOR ax('x',N_F);

        for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
        {
            ax[j_pivotz-N_C] = G_tau(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
        }

        ADDEND_BETH(ax);
    }

    FN_EXIT_POINT;
}


/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/


L_DOUBLE SVdata::opt_actmedium_get_e(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );

    L_DOUBLE result;

    result = fast_e_tau[iz];

    if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] == 2 ) )
    {
        result -= fast_rho[iz];
    }

    else if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] == 1 ) )
    {
        result += fast_rho_star[iz];
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_actmedium_get_e_pos(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] >= 0 );

    L_DOUBLE result;

    result = fast_e_tau[iz];

    if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] != 2 ) )
    {
        result += fast_rho[iz];
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_actmedium_get_e_neg(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] <= 0 );

    L_DOUBLE result;

    result = fast_e_tau[iz];

    if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] != 1 ) )
    {
        result -= fast_rho_star[iz];
    }

    FN_EXIT_POINT result;
}


L_DOUBLE SVdata::opt_actmedium_get_G_tau(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    L_DOUBLE result;
    long i_pivotz,j_pivotz;

    i_pivotz = fast_minv[iz]-1;
    j_pivotz = fast_minv[jz]-1;

    if ( ( i_pivotz >= N_C ) && ( j_pivotz >= N_C ) )
    {
        result = G_tau(i_pivotz-N_C,j_pivotz-N_C);
    }

    else
    {
        result = opt_actmedium_calc_G_tau_from_scratch(iz,jz);
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_actmedium_get_G_extras(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    FN_EXIT_POINT 0.0;

    iz = 1;
    jz = 1;
}

L_DOUBLE SVdata::opt_actmedium_get_G_tau_pivot(long i_pivotz, long j_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= 0 );
    THROW_ASSERT( i_pivotz <  N );
    THROW_ASSERT( j_pivotz >= 0 );
    THROW_ASSERT( j_pivotz <  N );

    L_DOUBLE result;

    if ( ( i_pivotz >= N_C ) && ( j_pivotz >= N_C ) )
    {
        result = G_tau(i_pivotz-N_C,j_pivotz-N_C);
    }

    else
    {
        result = opt_actmedium_calc_G_tau_from_scratch(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_actmedium_calc_G_tau_from_scratch(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    L_DOUBLE result;

    result = K.kernel(*(fast_x[iz]),*(fast_x[jz]),fast_opt_x[iz],fast_opt_x[jz],0.0,0.0,iz+1,jz+1,fast_z[iz],fast_z[jz]);

    FN_EXIT_POINT result;
}


void SVdata::opt_actmedium_step_none(L_DOUBLE d_b)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    d_b = 0.0;

    FN_EXIT_POINT;
}

void SVdata::opt_actmedium_step_one(L_DOUBLE d_alpha, long iz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    d_alpha = 0.0;
    iz      = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_actmedium_step_two(L_DOUBLE d_b, L_DOUBLE d_f, L_DOUBLE d_alpha1, L_DOUBLE d_alpha2, long i1z, long i2z)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    d_b      = 0.0;
    d_f      = 0.0;
    d_alpha1 = 0.0;
    d_alpha2 = 0.0;
    i1z      = 0;
    i2z      = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_actmedium_step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( d_alpha_pivot.get_effective_size() == d_e_pivot.get_effective_size() );
    THROW_ASSERT( d_alpha_pivot.get_effective_size() == finish_pivot-start_pivot+1 );
    THROW_ASSERT( start_pivot >= N_C+1 );
    THROW_ASSERT( start_pivot <= N );
    THROW_ASSERT( finish_pivot >= N_C+1 );
    THROW_ASSERT( finish_pivot <= N );
    THROW_ASSERT( start_pivot <= finish_pivot );

    long i_pivot,j_pivot;

    if ( n_ignore > 0 )
    {
        for ( i_pivot = start_pivot ; i_pivot < start_pivot+n_ignore ; i_pivot++ )
        {
            fast_alpha[fast_m[i_pivot-1]-1] += d_alpha_pivot[i_pivot-start_pivot];
        }
    }

    if ( start_pivot+n_ignore <= finish_pivot )
    {
        for ( i_pivot = start_pivot+n_ignore ; i_pivot <= finish_pivot ; i_pivot++ )
        {
            fast_alpha[fast_m[i_pivot-1]-1] += d_alpha_pivot[i_pivot-start_pivot];
            fast_e_tau[fast_m[i_pivot-1]-1] += d_e_pivot[i_pivot-start_pivot];
        }
    }

    f += d_f;

    if ( start_pivot > 1 )
    {
        for ( i_pivot = 1 ; i_pivot < start_pivot ; i_pivot++ )
        {
            for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
            {
                fast_e_tau[fast_m[i_pivot-1]-1] += opt_actmedium_get_G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
            }

            fast_e_tau[fast_m[i_pivot-1]-1] += d_b;
        }
    }

    if ( finish_pivot < N )
    {
        for ( i_pivot = finish_pivot+1 ; i_pivot <= N ; i_pivot++ )
        {
            for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
            {
                fast_e_tau[fast_m[i_pivot-1]-1] += opt_actmedium_get_G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
            }

            fast_e_tau[fast_m[i_pivot-1]-1] += d_b;
        }
    }

    b += d_b;

    opt_fix_e_free_pivot();

    FN_EXIT_POINT;
}

void SVdata::opt_actmedium_step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( d_alpha_pivot.get_effective_size() == finish_pivot-start_pivot+1 );
    THROW_ASSERT( start_pivot >= N_C+1 );
    THROW_ASSERT( start_pivot <= N );
    THROW_ASSERT( finish_pivot >= N_C+1 );
    THROW_ASSERT( finish_pivot <= N );
    THROW_ASSERT( start_pivot <= finish_pivot );

    long i_pivot,j_pivot;

    for ( i_pivot = start_pivot ; i_pivot <= finish_pivot ; i_pivot++ )
    {
        fast_alpha[fast_m[i_pivot-1]-1] += d_alpha_pivot[i_pivot-start_pivot];
        f += d_alpha_pivot[i_pivot-start_pivot];
    }

    for ( i_pivot = 1 ; i_pivot <= N ; i_pivot++ )
    {
        for ( j_pivot = start_pivot ; j_pivot <= finish_pivot ; j_pivot++ )
        {
            fast_e_tau[fast_m[i_pivot-1]-1] += opt_actmedium_get_G_tau(fast_m[i_pivot-1]-1,fast_m[j_pivot-1]-1)*d_alpha_pivot[j_pivot-start_pivot];
        }

        fast_e_tau[fast_m[i_pivot-1]-1] += d_b;
    }

    b += d_b;

    opt_fix_e_free_pivot();

    FN_EXIT_POINT;
}

void SVdata::opt_actmedium_constrain_Z(long iz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    iz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_actmedium_constrain_L(long iz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    iz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_actmedium_constrain_U(long iz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    iz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_actmedium_free_L(long iz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    iz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_actmedium_free_U(long iz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    iz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_actmedium_constrain_Z_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= N_C );
    THROW_ASSERT( i_pivotz <  N   );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != -2 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != +2 );

    long iz,j_pivotz;

    if ( fast_tau[fast_m[i_pivotz]-1] != 0 )
    {
        e_free_pivot.remove(i_pivotz+1-N_C);

        iz = fast_m[i_pivotz]-1;

        SHRINK_BETH(i_pivotz+1-N_C);

        G_tau.remove(i_pivotz+1-N_C);

        if ( fast_tau[iz] == -1 )
        {
            if ( fast_contype[iz] != 1 )
            {
                fast_e_tau[iz] += fast_rho_star[iz];
            }

            N_FN--;
        }

        else
        {
            if ( fast_contype[iz] != 2 )
            {
                fast_e_tau[iz] -= fast_rho[iz];
            }

            N_FP--;
        }

        fast_tau[iz] = 0;

        m.bswap(N_Z+1,i_pivotz+1);

        for ( j_pivotz = N_Z ; j_pivotz <= i_pivotz ; j_pivotz++ )
        {
            fast_minv[fast_m[j_pivotz]-1] = j_pivotz+1;
        }

        N_Z++;
        N_F--;
        N_C++;
        N_S--;
    }

    FN_EXIT_POINT;
}

void SVdata::opt_actmedium_constrain_L_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= N_C );
    THROW_ASSERT( i_pivotz <  N   );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != 0 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != +1 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != +2 );

    long iz,j_pivotz;

    if ( fast_tau[fast_m[i_pivotz]-1] != -2 )
    {
        e_free_pivot.remove(i_pivotz+1-N_C);

        iz = fast_m[i_pivotz]-1;

        SHRINK_BETH(i_pivotz+1-N_C);

        G_tau.remove(i_pivotz+1-N_C);

        fast_tau[iz] = -2;

        m.bswap(N_Z+N_L+1,i_pivotz+1);

        for ( j_pivotz = N_Z+N_L ; j_pivotz <= i_pivotz ; j_pivotz++ )
        {
            fast_minv[fast_m[j_pivotz]-1] = j_pivotz+1;
        }

        N_L++;
        N_F--;
        N_FN--;
        N_C++;
    }

    FN_EXIT_POINT;
}

void SVdata::opt_actmedium_constrain_U_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= N_C );
    THROW_ASSERT( i_pivotz <  N   );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != -2 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != -1 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != 0 );

    long iz,j_pivotz;

    if ( fast_tau[fast_m[i_pivotz]-1] != +2 )
    {
        e_free_pivot.remove(i_pivotz+1-N_C);

        iz = fast_m[i_pivotz]-1;

        SHRINK_BETH(i_pivotz+1-N_C);

        G_tau.remove(i_pivotz+1-N_C);

        fast_tau[iz] = +2;

        m.bswap(N_C+1,i_pivotz+1);

        for ( j_pivotz = N_C ; j_pivotz <= i_pivotz ; j_pivotz++ )
        {
            fast_minv[fast_m[j_pivotz]-1] = j_pivotz+1;
        }

        N_U++;
        N_F--;
        N_FP--;
        N_C++;
    }

    FN_EXIT_POINT;
}

void SVdata::opt_actmedium_free_L_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i_pivotz >= 0 );
    THROW_ASSERT( i_pivotz <= N_Z+N_L );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != 1 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != 2 );
    THROW_ASSERT( fast_contype[fast_m[i_pivotz]-1] != 0 );
    THROW_ASSERT( fast_contype[fast_m[i_pivotz]-1] != 2 );

    long iz,j_pivotz;
    int old_tau;

    if ( fast_tau[fast_m[i_pivotz]-1] != -1 )
    {
        iz = fast_m[i_pivotz]-1;

        old_tau = fast_tau[iz];

        if ( old_tau == 0 )
        {
            N_Z--;
            N_F++;
            N_FN++;
            N_C--;
            N_S++;
        }

        else
        {
            N_L--;
            N_F++;
            N_FN++;
            N_C--;
        }

        fast_tau[iz] = -1;

        m.fswap(i_pivotz+1,N);

        for ( j_pivotz = i_pivotz ; j_pivotz < N ; j_pivotz++ )
        {
            fast_minv[fast_m[j_pivotz]-1] = j_pivotz+1;
        }

        i_pivotz = N-1;

        if ( old_tau == 0 )
        {
            fVECTOR ax('x',N_F);

            for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
            {
                ax[j_pivotz-N_C] = opt_actmedium_calc_G_tau_from_scratch(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
            }

            G_tau.addend(ax[N_F-1],ax);

            if ( fast_contype[iz] != 1 )
            {
                fast_e_tau[iz] -= fast_rho_star[iz];
            }
        }

        else
        {
            fVECTOR ax('x',N_F);

            for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
            {
                ax[j_pivotz-N_C] = opt_actmedium_calc_G_tau_from_scratch(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
            }

            G_tau.addend(ax[N_F-1],ax);
        }

        e_free_pivot.addend(fast_e_tau[iz]);

        fVECTOR ax('x',N_F);

        for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
        {
            ax[j_pivotz-N_C] = G_tau(i_pivotz-N_C,j_pivotz-N_C);
        }

        ADDEND_BETH(ax);
    }

    FN_EXIT_POINT;
}

void SVdata::opt_actmedium_free_U_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( ( ( i_pivotz >= 0 ) && ( i_pivotz < N_Z ) ) || ( ( i_pivotz >= N_Z+N_L ) && ( i_pivotz < N_C ) ) );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != -1 );
    THROW_ASSERT( fast_tau[fast_m[i_pivotz]-1] != -2 );
    THROW_ASSERT( fast_contype[fast_m[i_pivotz]-1] != 0 );
    THROW_ASSERT( fast_contype[fast_m[i_pivotz]-1] != 1 );

    long iz,j_pivotz;
    int old_tau;

    if ( fast_tau[fast_m[i_pivotz]-1] != +1 )
    {
        iz = fast_m[i_pivotz]-1;

        old_tau = fast_tau[iz];

        if ( old_tau == 0 )
        {
            N_Z--;
            N_F++;
            N_FP++;
            N_C--;
            N_S++;
        }

        else
        {
            N_U--;
            N_F++;
            N_FP++;
            N_C--;
        }

        fast_tau[iz] = +1;

        m.fswap(i_pivotz+1,N);

        for ( j_pivotz = i_pivotz ; j_pivotz < N ; j_pivotz++ )
        {
            fast_minv[fast_m[j_pivotz]-1] = j_pivotz+1;
        }

        i_pivotz = N-1;

        if ( old_tau == 0 )
        {
            fVECTOR ax('x',N_F);

            for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
            {
                ax[j_pivotz-N_C] = opt_actmedium_calc_G_tau_from_scratch(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
            }

            G_tau.addend(ax[N_F-1],ax);

            if ( fast_contype[iz] != 2 )
            {
                fast_e_tau[iz] += fast_rho[iz];
            }
        }

        else
        {
            fVECTOR ax('x',N_F);

            for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
            {
                ax[j_pivotz-N_C] = opt_actmedium_calc_G_tau_from_scratch(fast_m[i_pivotz]-1,fast_m[j_pivotz]-1);
            }

            G_tau.addend(ax[N_F-1],ax);
        }

        e_free_pivot.addend(fast_e_tau[iz]);

        fVECTOR ax('x',N_F);

        for ( j_pivotz = N_C ; j_pivotz < N ; j_pivotz++ )
        {
            ax[j_pivotz-N_C] = G_tau(i_pivotz-N_C,j_pivotz-N_C);
        }

        ADDEND_BETH(ax);
    }

    FN_EXIT_POINT;
}


/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/


#ifdef ROWPLUS
#undef ROWPLUS
#endif
#define ROWPLUS(_isrow_,_iz_,_jz_) ((_isrow_)[(_iz_)])

L_DOUBLE SVdata::opt_d2csmo_get_e(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );

    L_DOUBLE result;

    result = fast_e_tau[iz];

    if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] == 2 ) )
    {
        result -= fast_rho[iz];
    }

    else if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] == 1 ) )
    {
        result += fast_rho_star[iz];
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_d2csmo_get_e_pos(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] >= 0 );

    L_DOUBLE result;

    result = fast_e_tau[iz];

    if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] != 2 ) )
    {
        result += fast_rho[iz];
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_d2csmo_get_e_neg(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] <= 0 );

    L_DOUBLE result;

    result = fast_e_tau[iz];

    if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] != 1 ) )
    {
        result -= fast_rho_star[iz];
    }

    FN_EXIT_POINT result;
}


L_DOUBLE SVdata::opt_d2csmo_get_G_tau(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    L_DOUBLE result;

    result = (kerncache).opt_getval_statcache(this,iz+1,jz+1);

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_d2csmo_get_G_extras(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    FN_EXIT_POINT 0.0;

    iz = 1;
    jz = 1;
}

L_DOUBLE SVdata::opt_d2csmo_get_G_tau_pivot(long i_pivotz, long j_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;
    j_pivotz = 0;

    FN_EXIT_POINT 0.0;
}

L_DOUBLE SVdata::opt_d2csmo_calc_G_tau_from_scratch(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    L_DOUBLE result;

    result = K.kernel(*(fast_x[iz]),*(fast_x[jz]),fast_opt_x[iz],fast_opt_x[jz],0.0,0.0,iz+1,jz+1,fast_z[iz],fast_z[jz]);

    FN_EXIT_POINT result;
}


void SVdata::opt_d2csmo_step_none(L_DOUBLE d_b)
{
    FN_ENTRY_POINT;

    long jz;

    for ( jz = 0 ; jz < N ; jz++ )
    {
        fast_e_tau[jz] += d_b;
    }

    b += d_b;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_step_one(L_DOUBLE d_alpha, long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( ( fast_tau[iz] == -1 ) || ( fast_tau[iz] == +1 ) );

    long jz;
    L_DOUBLE *isrow;

    fast_alpha[iz] += d_alpha;

    f += d_alpha;

    isrow = KERNGETROWOPT(iz);

    for ( jz = 0 ; jz < N ; jz++ )
    {
        fast_e_tau[jz] += ROWPLUS(isrow,jz,iz) * d_alpha;
    }

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_step_two(L_DOUBLE d_b, L_DOUBLE d_f, L_DOUBLE d_alpha1, L_DOUBLE d_alpha2, long i1z, long i2z)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i1z >= 0 );
    THROW_ASSERT( i1z <  N );
    THROW_ASSERT( i2z >= 0 );
    THROW_ASSERT( i2z <  N );
    THROW_ASSERT( ( fast_tau[i1z] == -1 ) || ( fast_tau[i1z] == +1 ) );
    THROW_ASSERT( ( fast_tau[i2z] == -1 ) || ( fast_tau[i2z] == +1 ) );

    long jz;
    L_DOUBLE *isrow1;
    L_DOUBLE *isrow2;

    fast_alpha[i1z] += d_alpha1;
    fast_alpha[i2z] += d_alpha2;

    f += d_f;

    isrow1 = KERNGETROWOPT(i1z);
    isrow2 = KERNGETROWOPT(i2z);

    for ( jz = 0 ; jz < N ; jz++ )
    {
        fast_e_tau[jz] += ROWPLUS(isrow1,jz,i1z) * d_alpha1;
        fast_e_tau[jz] += ROWPLUS(isrow2,jz,i2z) * d_alpha2;
        fast_e_tau[jz] += d_b;
    }

    b += d_b;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    start_pivot   = 0;
    finish_pivot  = 0;
    d_alpha_pivot = 0.0;
    d_b           = 0.0;
    d_e_pivot     = 0.0;
    d_f           = 0.0;
    n_ignore      = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    start_pivot   = 0;
    finish_pivot  = 0;
    d_alpha_pivot = 0.0;
    d_b           = 0.0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_constrain_Z(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != -2 );
    THROW_ASSERT( fast_tau[iz] != +2 );

    if ( ( fast_tau[iz] == -1 ) && ( fast_contype[iz] != 1 ) )
    {
        fast_e_tau[iz] += fast_rho_star[iz];
    }

    else if ( ( fast_tau[iz] == +1 ) && ( fast_contype[iz] != 2 ) )
    {
        fast_e_tau[iz] -= fast_rho[iz];
    }

    fast_tau[iz] = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_constrain_L(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != 0 );
    THROW_ASSERT( fast_tau[iz] != +1 );
    THROW_ASSERT( fast_tau[iz] != +2 );

    fast_tau[iz] = -2;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_constrain_U(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != 0 );
    THROW_ASSERT( fast_tau[iz] != -1 );
    THROW_ASSERT( fast_tau[iz] != -2 );

    fast_tau[iz] = +2;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_free_L(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != +1 );
    THROW_ASSERT( fast_tau[iz] != +2 );
    THROW_ASSERT( fast_contype[iz] !=  0 );
    THROW_ASSERT( fast_contype[iz] != +2 );

    if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] != 1 ) )
    {
        fast_e_tau[iz] -= fast_rho_star[iz];
    }

    fast_tau[iz] = -1;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_free_U(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != -1 );
    THROW_ASSERT( fast_tau[iz] != -2 );
    THROW_ASSERT( fast_contype[iz] !=  0 );
    THROW_ASSERT( fast_contype[iz] != +1 );

    if ( ( fast_tau[iz] == 0 ) && ( fast_contype[iz] != 2 ) )
    {
        fast_e_tau[iz] += fast_rho[iz];
    }

    fast_tau[iz] = +1;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_constrain_Z_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_constrain_L_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_constrain_U_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_free_L_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_free_U_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;

    FN_EXIT_POINT;
}


/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/


#ifdef ROWPLUS
#undef ROWPLUS
#endif
#define ROWPLUS(_isrow_,_iz_,_jz_) ((_isrow_)[(_iz_)])

L_DOUBLE SVdata::opt_d2csmo_pattern_get_e(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );

    L_DOUBLE result;

    result = fast_e_tau[iz];

    if ( fast_tau[iz] == 0 )
    {
        if ( fast_contype[iz] == 2 )
        {
            result -= fast_rho[iz];
        }

        else
        {
            result += fast_rho_star[iz];
        }
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_d2csmo_pattern_get_e_pos(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] >= 0 );

    FN_EXIT_POINT fast_e_tau[iz];
}

L_DOUBLE SVdata::opt_d2csmo_pattern_get_e_neg(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] <= 0 );

    FN_EXIT_POINT fast_e_tau[iz];
}


L_DOUBLE SVdata::opt_d2csmo_pattern_get_G_tau(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    L_DOUBLE result;

    result = (kerncache).opt_getval_statcache(this,iz+1,jz+1);

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_d2csmo_pattern_get_G_extras(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    FN_EXIT_POINT 0.0;

    iz = 1;
    jz = 1;
}

L_DOUBLE SVdata::opt_d2csmo_pattern_get_G_tau_pivot(long i_pivotz, long j_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;
    j_pivotz = 0;

    FN_EXIT_POINT 0.0;
}

L_DOUBLE SVdata::opt_d2csmo_pattern_calc_G_tau_from_scratch(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    L_DOUBLE result;

    result = K.kernel(*(fast_x[iz]),*(fast_x[jz]),fast_opt_x[iz],fast_opt_x[jz],0.0,0.0,iz+1,jz+1,fast_z[iz],fast_z[jz]);

    FN_EXIT_POINT result;
}


void SVdata::opt_d2csmo_pattern_step_none(L_DOUBLE d_b)
{
    FN_ENTRY_POINT;

    long jz;

    for ( jz = 0 ; jz < N ; jz++ )
    {
        fast_e_tau[jz] += d_b;
    }

    b += d_b;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_step_one(L_DOUBLE d_alpha, long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( ( fast_tau[iz] == -1 ) || ( fast_tau[iz] == +1 ) );

    long jz;
    L_DOUBLE *isrow;

    fast_alpha[iz] += d_alpha;

    f += d_alpha;

    isrow = KERNGETROWOPT(iz);

    for ( jz = 0 ; jz < N ; jz++ )
    {
        fast_e_tau[jz] += ROWPLUS(isrow,jz,iz) * d_alpha;
    }

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_step_two(L_DOUBLE d_b, L_DOUBLE d_f, L_DOUBLE d_alpha1, L_DOUBLE d_alpha2, long i1z, long i2z)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i1z >= 0 );
    THROW_ASSERT( i1z <  N );
    THROW_ASSERT( i2z >= 0 );
    THROW_ASSERT( i2z <  N );
    THROW_ASSERT( ( fast_tau[i1z] == -1 ) || ( fast_tau[i1z] == +1 ) );
    THROW_ASSERT( ( fast_tau[i2z] == -1 ) || ( fast_tau[i2z] == +1 ) );

    long jz;
    L_DOUBLE *isrow1;
    L_DOUBLE *isrow2;

    fast_alpha[i1z] += d_alpha1;
    fast_alpha[i2z] += d_alpha2;

    f += d_f;

    isrow1 = KERNGETROWOPT(i1z);
    isrow2 = KERNGETROWOPT(i2z);

    for ( jz = 0 ; jz < N ; jz++ )
    {
        fast_e_tau[jz] += ROWPLUS(isrow1,jz,i1z) * d_alpha1;
        fast_e_tau[jz] += ROWPLUS(isrow2,jz,i2z) * d_alpha2;
        fast_e_tau[jz] += d_b;
    }

    b += d_b;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    start_pivot   = 0;
    finish_pivot  = 0;
    d_alpha_pivot = 0.0;
    d_b           = 0.0;
    d_e_pivot     = 0.0;
    d_f           = 0.0;
    n_ignore      = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    start_pivot   = 0;
    finish_pivot  = 0;
    d_alpha_pivot = 0.0;
    d_b           = 0.0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_constrain_Z(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != -2 );
    THROW_ASSERT( fast_tau[iz] != +2 );

    fast_tau[iz] = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_constrain_L(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != 0 );
    THROW_ASSERT( fast_tau[iz] != +1 );
    THROW_ASSERT( fast_tau[iz] != +2 );

    fast_tau[iz] = -2;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_constrain_U(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != 0 );
    THROW_ASSERT( fast_tau[iz] != -1 );
    THROW_ASSERT( fast_tau[iz] != -2 );

    fast_tau[iz] = +2;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_free_L(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != +1 );
    THROW_ASSERT( fast_tau[iz] != +2 );
    THROW_ASSERT( fast_contype[iz] !=  0 );
    THROW_ASSERT( fast_contype[iz] != +2 );

    fast_tau[iz] = -1;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_free_U(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != -1 );
    THROW_ASSERT( fast_tau[iz] != -2 );
    THROW_ASSERT( fast_contype[iz] !=  0 );
    THROW_ASSERT( fast_contype[iz] != +1 );

    fast_tau[iz] = +1;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_constrain_Z_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_constrain_L_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_constrain_U_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_free_L_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_free_U_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;

    FN_EXIT_POINT;
}



















/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/


#ifdef ROWPLUS
#undef ROWPLUS
#endif
#define ROWPLUS(_isrow_,_iz_,_jz_) ((_isrow_)[(_iz_)])

#undef KERNGETROWOPT
#define KERNGETROWOPT(__i) (kerncache.opt_getrow_fast(this,((__i)+1)))

L_DOUBLE SVdata::opt_d2csmo_pattern_kfast_get_e(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );

    L_DOUBLE result;

    result = fast_e_tau[iz];

    if ( fast_tau[iz] == 0 )
    {
        if ( fast_contype[iz] == 2 )
        {
            result -= fast_rho[iz];
        }

        else
        {
            result += fast_rho_star[iz];
        }
    }

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_d2csmo_pattern_kfast_get_e_pos(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] >= 0 );

    FN_EXIT_POINT fast_e_tau[iz];
}

L_DOUBLE SVdata::opt_d2csmo_pattern_kfast_get_e_neg(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] <= 0 );

    FN_EXIT_POINT fast_e_tau[iz];
}


L_DOUBLE SVdata::opt_d2csmo_pattern_kfast_get_G_tau(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    L_DOUBLE result;

    result = fast_kernel_fn(fast_opt_x[iz],fast_opt_x[jz],fast_uc,fast_uv,(double *) fast_covw);

    FN_EXIT_POINT result;
}

L_DOUBLE SVdata::opt_d2csmo_pattern_kfast_get_G_extras(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    FN_EXIT_POINT 0.0;

    iz = 1;
    jz = 1;
}

L_DOUBLE SVdata::opt_d2csmo_pattern_kfast_get_G_tau_pivot(long i_pivotz, long j_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;
    j_pivotz = 0;

    FN_EXIT_POINT 0.0;
}

L_DOUBLE SVdata::opt_d2csmo_pattern_kfast_calc_G_tau_from_scratch(long iz, long jz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( jz >= 0 );
    THROW_ASSERT( jz <  N );

    L_DOUBLE result;

    result = K.kernel(*(fast_x[iz]),*(fast_x[jz]),fast_opt_x[iz],fast_opt_x[jz],0.0,0.0,iz+1,jz+1,fast_z[iz],fast_z[jz]);

    FN_EXIT_POINT result;
}


void SVdata::opt_d2csmo_pattern_kfast_step_none(L_DOUBLE d_b)
{
    FN_ENTRY_POINT;

    long jz;

    for ( jz = 0 ; jz < N ; jz++ )
    {
        fast_e_tau[jz] += d_b;
    }

    b += d_b;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_kfast_step_one(L_DOUBLE d_alpha, long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( ( fast_tau[iz] == -1 ) || ( fast_tau[iz] == +1 ) );

    long jz;
    L_DOUBLE *isrow;

    fast_alpha[iz] += d_alpha;

    f += d_alpha;

    isrow = KERNGETROWOPT(iz);

    for ( jz = 0 ; jz < N ; jz++ )
    {
        fast_e_tau[jz] += ROWPLUS(isrow,jz,iz) * d_alpha;
    }

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_kfast_step_two(L_DOUBLE d_b, L_DOUBLE d_f, L_DOUBLE d_alpha1, L_DOUBLE d_alpha2, long i1z, long i2z)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( i1z >= 0 );
    THROW_ASSERT( i1z <  N );
    THROW_ASSERT( i2z >= 0 );
    THROW_ASSERT( i2z <  N );
    THROW_ASSERT( ( fast_tau[i1z] == -1 ) || ( fast_tau[i1z] == +1 ) );
    THROW_ASSERT( ( fast_tau[i2z] == -1 ) || ( fast_tau[i2z] == +1 ) );

    long jz;
    L_DOUBLE *isrow1;
    L_DOUBLE *isrow2;

    fast_alpha[i1z] += d_alpha1;
    fast_alpha[i2z] += d_alpha2;

    f += d_f;

    isrow1 = KERNGETROWOPT(i1z);
    isrow2 = KERNGETROWOPT(i2z);

    for ( jz = 0 ; jz < N ; jz++ )
    {
        fast_e_tau[jz] += ROWPLUS(isrow1,jz,i1z) * d_alpha1;
        fast_e_tau[jz] += ROWPLUS(isrow2,jz,i2z) * d_alpha2;
        fast_e_tau[jz] += d_b;
    }

    b += d_b;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_kfast_step_general_pivot(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b, fVECTOR &d_e_pivot, L_DOUBLE d_f, long n_ignore)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    start_pivot   = 0;
    finish_pivot  = 0;
    d_alpha_pivot = 0.0;
    d_b           = 0.0;
    d_e_pivot     = 0.0;
    d_f           = 0.0;
    n_ignore      = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_kfast_step_general_pivotx(long start_pivot, long finish_pivot, fVECTOR &d_alpha_pivot, L_DOUBLE d_b)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    start_pivot   = 0;
    finish_pivot  = 0;
    d_alpha_pivot = 0.0;
    d_b           = 0.0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_kfast_constrain_Z(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != -2 );
    THROW_ASSERT( fast_tau[iz] != +2 );

    fast_tau[iz] = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_kfast_constrain_L(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != 0 );
    THROW_ASSERT( fast_tau[iz] != +1 );
    THROW_ASSERT( fast_tau[iz] != +2 );

    fast_tau[iz] = -2;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_kfast_constrain_U(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != 0 );
    THROW_ASSERT( fast_tau[iz] != -1 );
    THROW_ASSERT( fast_tau[iz] != -2 );

    fast_tau[iz] = +2;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_kfast_free_L(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != +1 );
    THROW_ASSERT( fast_tau[iz] != +2 );
    THROW_ASSERT( fast_contype[iz] !=  0 );
    THROW_ASSERT( fast_contype[iz] != +2 );

    fast_tau[iz] = -1;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_kfast_free_U(long iz)
{
    FN_ENTRY_POINT;

    THROW_ASSERT( iz >= 0 );
    THROW_ASSERT( iz <  N );
    THROW_ASSERT( fast_tau[iz] != -1 );
    THROW_ASSERT( fast_tau[iz] != -2 );
    THROW_ASSERT( fast_contype[iz] !=  0 );
    THROW_ASSERT( fast_contype[iz] != +1 );

    fast_tau[iz] = +1;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_kfast_constrain_Z_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_kfast_constrain_L_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_kfast_constrain_U_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_kfast_free_L_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;

    FN_EXIT_POINT;
}

void SVdata::opt_d2csmo_pattern_kfast_free_U_pivot(long i_pivotz)
{
    FN_ENTRY_POINT;

    L_THROW(0);

    i_pivotz = 0;

    FN_EXIT_POINT;
}



