
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


/*
   Note that the caching, factorisation and optimisation options cannot be
   combined - that is, they are exclusive.  The "hint" flags, however, can
   be combined as required to form a composite set.
*/

#define E_CACHE_FLAG_MASK               0x00000000f
#define E_CACHE_FULL_FLAG               0x000000001
#define E_CACHE_SUPP_FLAG               0x000000002
#define E_CACHE_FREE_FLAG               0x000000003
/* fallback no e cache */

#define KERN_CACHE_FLAG_MASK            0x0000000f0
#define KERN_CACHE_FULL_FLAG            0x000000010
#define KERN_CACHE_FREE_FLAG            0x000000020
#define KERN_CACHE_DYNA_FLAG            0x000000030
/* fallback no kernel cache */

#define H_FACT_FLAG_MASK                0x000000f00
#define H_FACT_INVE_FLAG                0x000000100
#define H_FACT_CHOL_FLAG                0x000000200
#define H_FACT_SLOW_FLAG                0x000000300
/* fallback no factorisation */

#define OPT_TYPE_FLAG_MASK              0x00000f000
#define OPT_ACTIVE_SET_FLAG             0x000001000
#define OPT_PLATT_SMO_FLAG              0x000002000
#define OPT_DANIEL_D2C_FLAG             0x000004000
/* fallback platts smo method */

#define SVM_HINT_MASK                   0x000ff0000
#define SVM_HINT_RHO_ZERO_FLAG          0x000010000
#define SVM_HINT_GAMMA_ZERO_FLAG        0x000020000
#define SVM_HINT_MU_ZERO_FLAG           0x000040000
#define SVM_HINT_PATTERN_FLAG           0x000080000

#define DEFAULT_SVFLAGS                 0x000070000


#define CLEAR_FLAGS(svflags)                    { (svflags) = DEFAULT_SVFLAGS; }

#define MAKE_E_CACHE_NONE(svflags)              { (svflags) =   (svflags) & (~E_CACHE_FLAG_MASK); }
#define MAKE_E_CACHE_FULL(svflags)              { (svflags) = ( (svflags) & (~E_CACHE_FLAG_MASK) ) | E_CACHE_FULL_FLAG; }
#define MAKE_E_CACHE_SUPP(svflags)              { (svflags) = ( (svflags) & (~E_CACHE_FLAG_MASK) ) | E_CACHE_SUPP_FLAG; }
#define MAKE_E_CACHE_FREE(svflags)              { (svflags) = ( (svflags) & (~E_CACHE_FLAG_MASK) ) | E_CACHE_FREE_FLAG; }

#define MAKE_KERN_CACHE_NONE(svflags)           { (svflags) =   (svflags) & (~KERN_CACHE_FLAG_MASK); }
#define MAKE_KERN_CACHE_FULL(svflags)           { (svflags) = ( (svflags) & (~KERN_CACHE_FLAG_MASK) ) | KERN_CACHE_FULL_FLAG; }
#define MAKE_KERN_CACHE_FREE(svflags)           { (svflags) = ( (svflags) & (~KERN_CACHE_FLAG_MASK) ) | KERN_CACHE_FREE_FLAG; }
#define MAKE_KERN_CACHE_DYNA(svflags)           { (svflags) = ( (svflags) & (~KERN_CACHE_FLAG_MASK) ) | KERN_CACHE_DYNA_FLAG; }

#define MAKE_H_FACT_NONE(svflags)               { (svflags) =   (svflags) & (~H_FACT_FLAG_MASK); }
#define MAKE_H_FACT_INVE(svflags)               { (svflags) = ( (svflags) & (~H_FACT_FLAG_MASK) ) | H_FACT_INVE_FLAG; }
#define MAKE_H_FACT_CHOL(svflags)               { (svflags) = ( (svflags) & (~H_FACT_FLAG_MASK) ) | H_FACT_CHOL_FLAG; }
#define MAKE_H_FACT_SLOW(svflags)               { (svflags) = ( (svflags) & (~H_FACT_FLAG_MASK) ) | H_FACT_SLOW_FLAG; }

#define MAKE_OPT_NONE(svflags)                  { (svflags) =   (svflags) & (~OPT_TYPE_FLAG_MASK); }
#define MAKE_OPT_ACTIVE_SET(svflags)            { (svflags) = ( (svflags) & (~OPT_TYPE_FLAG_MASK) ) | OPT_ACTIVE_SET_FLAG; }
#define MAKE_OPT_PLATT_SMO(svflags)             { (svflags) = ( (svflags) & (~OPT_TYPE_FLAG_MASK) ) | OPT_PLATT_SMO_FLAG;  }
#define MAKE_OPT_DANIEL_D2C(svflags)            { (svflags) = ( (svflags) & (~OPT_TYPE_FLAG_MASK) ) | OPT_DANIEL_D2C_FLAG; }

#define MAKE_SVM_HINT_EMPTY(svflags)            { (svflags) =   (svflags) & (~SVM_HINT_MASK); }
#define MAKE_SVM_HINT_RHO_ZERO(svflags)         { (svflags) =   (svflags) | SVM_HINT_RHO_ZERO_FLAG;   }
#define MAKE_SVM_HINT_GAMMA_ZERO(svflags)       { (svflags) =   (svflags) | SVM_HINT_GAMMA_ZERO_FLAG; }
#define MAKE_SVM_HINT_MU_ZERO(svflags)          { (svflags) =   (svflags) | SVM_HINT_MU_ZERO_FLAG;    }
#define MAKE_SVM_HINT_PATTERN(svflags)          { (svflags) =   (svflags) | SVM_HINT_PATTERN_FLAG;    }
#define MAKE_SVM_HINT_RHO_NONZERO(svflags)      { (svflags) =   (svflags) & (~SVM_HINT_RHO_ZERO_FLAG);   }
#define MAKE_SVM_HINT_GAMMA_NONZERO(svflags)    { (svflags) =   (svflags) & (~SVM_HINT_GAMMA_ZERO_FLAG); }
#define MAKE_SVM_HINT_MU_NONZERO(svflags)       { (svflags) =   (svflags) & (~SVM_HINT_MU_ZERO_FLAG);    }
#define MAKE_SVM_HINT_NONPATTERN(svflags)       { (svflags) =   (svflags) & (~SVM_HINT_PATTERN_FLAG);    }


#define E_CACHE_FULL(svflags)           ( ( (svflags) & E_CACHE_FLAG_MASK ) == E_CACHE_FULL_FLAG )
#define E_CACHE_SUPP(svflags)           ( ( (svflags) & E_CACHE_FLAG_MASK ) == E_CACHE_SUPP_FLAG )
#define E_CACHE_FREE(svflags)           ( ( (svflags) & E_CACHE_FLAG_MASK ) == E_CACHE_FREE_FLAG )

#define KERN_CACHE_FULL(svflags)        ( ( (svflags) & KERN_CACHE_FLAG_MASK ) == KERN_CACHE_FULL_FLAG )
#define KERN_CACHE_FREE(svflags)        ( ( (svflags) & KERN_CACHE_FLAG_MASK ) == KERN_CACHE_FREE_FLAG )
#define KERN_CACHE_DYNA(svflags)        ( ( (svflags) & KERN_CACHE_FLAG_MASK ) == KERN_CACHE_DYNA_FLAG )

#define H_FACT_INVE(svflags)            ( ( (svflags) & H_FACT_FLAG_MASK ) == H_FACT_INVE_FLAG )
#define H_FACT_CHOL(svflags)            ( ( (svflags) & H_FACT_FLAG_MASK ) == H_FACT_CHOL_FLAG )
#define H_FACT_SLOW(svflags)            ( ( (svflags) & H_FACT_FLAG_MASK ) == H_FACT_SLOW_FLAG )

#define OPT_ACTIVE_SET(svflags)         ( ( (svflags) & OPT_TYPE_FLAG_MASK ) == OPT_ACTIVE_SET_FLAG )
#define OPT_PLATT_SMO(svflags)          ( ( (svflags) & OPT_TYPE_FLAG_MASK ) == OPT_PLATT_SMO_FLAG  )
#define OPT_DANIEL_D2C(svflags)         ( ( (svflags) & OPT_TYPE_FLAG_MASK ) == OPT_DANIEL_D2C_FLAG )

#define SVM_HINT_RHO_ZERO(svflags)      ( ( (svflags) & SVM_HINT_MASK ) & SVM_HINT_RHO_ZERO_FLAG   )
#define SVM_HINT_GAMMA_ZERO(svflags)    ( ( (svflags) & SVM_HINT_MASK ) & SVM_HINT_GAMMA_ZERO_FLAG )
#define SVM_HINT_MU_ZERO(svflags)       ( ( (svflags) & SVM_HINT_MASK ) & SVM_HINT_MU_ZERO_FLAG    )
#define SVM_HINT_PATTERN(svflags)       ( ( (svflags) & SVM_HINT_MASK ) & SVM_HINT_PATTERN_FLAG    )
#define SVM_HINT_EFREE_USED(svflags)    ( OPT_ACTIVE_SET(svflags) )
#define SVM_HINT_USE_PIVOTS(svflags)    ( OPT_ACTIVE_SET(svflags) || KERN_CACHE_FREE(svflags) )



/*
SMO/D2C optimised:            E_CACHE_FULL(svflags)
                              KERN_CACHE_DYNA(svflags)
                              OPT_DANIEL_D2C(svflags)
                              !SVM_HINT_RHO_ZERO(svflags)
                              SVM_HINT_GAMMA_ZERO(svflags)
                              SVM_HINT_MU_ZERO(svflags)
                              !SVM_HINT_EFREE_USED(svflags)
                              !SVM_HINT_USE_PIVOTS(svflags)

act small dataset optimised:  E_CACHE_FULL(svflags)
                              KERN_CACHE_FULL(svflags)
                              OPT_ACTIVE_SET(svflags)
                              !SVM_HINT_RHO_ZERO(svflags)
                              SVM_HINT_GAMMA_ZERO(svflags)
                              SVM_HINT_MU_ZERO(svflags)
                              SVM_HINT_EFREE_USED(svflags)
                              SVM_HINT_USE_PIVOTS(svflags)
                              some factorisation used

act medium dataset optimised: E_CACHE_FULL(svflags)
                              KERN_CACHE_FREE(svflags)
                              OPT_ACTIVE_SET(svflags)
                              !SVM_HINT_RHO_ZERO(svflags)
                              SVM_HINT_GAMMA_ZERO(svflags)
                              SVM_HINT_MU_ZERO(svflags)
                              SVM_HINT_EFREE_USED(svflags)
                              SVM_HINT_USE_PIVOTS(svflags)
                              some factorisation used
*/

#define D2CSMO_OPTIM(svflags)           ( E_CACHE_FULL(svflags) && KERN_CACHE_DYNA(svflags) && OPT_DANIEL_D2C(svflags) && !SVM_HINT_RHO_ZERO(svflags) && SVM_HINT_GAMMA_ZERO(svflags) && SVM_HINT_MU_ZERO(svflags) && !SVM_HINT_EFREE_USED(svflags) && !SVM_HINT_USE_PIVOTS(svflags) )
#define ACTIVE_SMALL_OPTIM(svflags)     ( E_CACHE_FULL(svflags) && KERN_CACHE_FULL(svflags) && OPT_ACTIVE_SET(svflags) && !SVM_HINT_RHO_ZERO(svflags) && SVM_HINT_GAMMA_ZERO(svflags) && SVM_HINT_MU_ZERO(svflags) &&  SVM_HINT_EFREE_USED(svflags) &&  SVM_HINT_USE_PIVOTS(svflags) && ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) ) )
#define ACTIVE_MEDIUM_OPTIM(svflags)    ( E_CACHE_FULL(svflags) && KERN_CACHE_FREE(svflags) && OPT_ACTIVE_SET(svflags) && !SVM_HINT_RHO_ZERO(svflags) && SVM_HINT_GAMMA_ZERO(svflags) && SVM_HINT_MU_ZERO(svflags) &&  SVM_HINT_EFREE_USED(svflags) &&  SVM_HINT_USE_PIVOTS(svflags) && ( H_FACT_INVE(svflags) || H_FACT_CHOL(svflags) ) )

/*
   Optimisation methods:

   These define the method of optimisation used by the SVM.  The basic
   format is SVM_OPT_www_xxx_yyy{_zzz}, where:

   www = optimisation method used.  Options are:

         ACTIVEST = active set method.
         PLATTSMO = Platt's SMO method.
         D2CALGOR = Daniel's D2C method.

   xxx = e (gradient) cache used.  Options are:

         FULLECACHE = cache all gradient values.
         SUPPECACHE = cache gradient values for support vectors.
         FREEECACHE = cache gradient values for unconstrained vectors.
         EMPTECACHE = no gradients cached.

   yyy = G (kernel) cached used.  Options are:

         FULLGCACHE = cache all kernel (hessian) matrix.
         FREEGCACHE = cache the unconstrained corner of the kernel matrix.
         DYNAGCACHE = dynamically cache the kernel with memory bounds.
                      (not usable when www = ACTIVE).
         EMPTGCACHE = no kernel values cached.
                      (not usable when www = ACTIVE).

   zzz = G factorisation used (relevant when www = ACTIVEST).  Options are:

         INVFACT  = inverse factorisation maintained.
         CHOLFACT = cholesky factorisation maintained.
         SLOWFACT = do all matrix inversions explicitly (very slow).
         NOFACT   = no factorisation kept.
*/

#define SVM_OPT_ACTIVEST_EMPTECACHE_FREEGCACHE_INVFACT  113
#define SVM_OPT_ACTIVEST_EMPTECACHE_FULLGCACHE_INVFACT  114
#define SVM_OPT_ACTIVEST_FREEECACHE_FREEGCACHE_INVFACT  123
#define SVM_OPT_ACTIVEST_FREEECACHE_FULLGCACHE_INVFACT  124
#define SVM_OPT_ACTIVEST_SUPPECACHE_FREEGCACHE_INVFACT  133
#define SVM_OPT_ACTIVEST_SUPPECACHE_FULLGCACHE_INVFACT  134
#define SVM_OPT_ACTIVEST_FULLECACHE_FREEGCACHE_INVFACT  143
#define SVM_OPT_ACTIVEST_FULLECACHE_FULLGCACHE_INVFACT  144
#define SVM_OPT_ACTIVEST_EMPTECACHE_FREEGCACHE_CHOLFACT 213
#define SVM_OPT_ACTIVEST_EMPTECACHE_FULLGCACHE_CHOLFACT 214
#define SVM_OPT_ACTIVEST_FREEECACHE_FREEGCACHE_CHOLFACT 223
#define SVM_OPT_ACTIVEST_FREEECACHE_FULLGCACHE_CHOLFACT 224
#define SVM_OPT_ACTIVEST_SUPPECACHE_FREEGCACHE_CHOLFACT 233
#define SVM_OPT_ACTIVEST_SUPPECACHE_FULLGCACHE_CHOLFACT 234
#define SVM_OPT_ACTIVEST_FULLECACHE_FREEGCACHE_CHOLFACT 243
#define SVM_OPT_ACTIVEST_FULLECACHE_FULLGCACHE_CHOLFACT 244
#define SVM_OPT_PLATTSMO_EMPTECACHE_EMPTGCACHE          311
#define SVM_OPT_PLATTSMO_EMPTECACHE_DYNAGCACHE          312
#define SVM_OPT_PLATTSMO_EMPTECACHE_FREEGCACHE          313
#define SVM_OPT_PLATTSMO_EMPTECACHE_FULLGCACHE          314
#define SVM_OPT_PLATTSMO_FREEECACHE_EMPTGCACHE          321
#define SVM_OPT_PLATTSMO_FREEECACHE_DYNAGCACHE          322
#define SVM_OPT_PLATTSMO_FREEECACHE_FREEGCACHE          323
#define SVM_OPT_PLATTSMO_FREEECACHE_FULLGCACHE          324
#define SVM_OPT_PLATTSMO_SUPPECACHE_EMPTGCACHE          331
#define SVM_OPT_PLATTSMO_SUPPECACHE_DYNAGCACHE          332
#define SVM_OPT_PLATTSMO_SUPPECACHE_FREEGCACHE          333
#define SVM_OPT_PLATTSMO_SUPPECACHE_FULLGCACHE          334
#define SVM_OPT_PLATTSMO_FULLECACHE_EMPTGCACHE          341
#define SVM_OPT_PLATTSMO_FULLECACHE_DYNAGCACHE          342
#define SVM_OPT_PLATTSMO_FULLECACHE_FREEGCACHE          343
#define SVM_OPT_PLATTSMO_FULLECACHE_FULLGCACHE          344
#define SVM_OPT_D2CALGOR_EMPTECACHE_EMPTGCACHE          411
#define SVM_OPT_D2CALGOR_EMPTECACHE_DYNAGCACHE          412
#define SVM_OPT_D2CALGOR_EMPTECACHE_FREEGCACHE          413
#define SVM_OPT_D2CALGOR_EMPTECACHE_FULLGCACHE          414
#define SVM_OPT_D2CALGOR_FREEECACHE_EMPTGCACHE          421
#define SVM_OPT_D2CALGOR_FREEECACHE_DYNAGCACHE          422
#define SVM_OPT_D2CALGOR_FREEECACHE_FREEGCACHE          423
#define SVM_OPT_D2CALGOR_FREEECACHE_FULLGCACHE          424
#define SVM_OPT_D2CALGOR_SUPPECACHE_EMPTGCACHE          431
#define SVM_OPT_D2CALGOR_SUPPECACHE_DYNAGCACHE          432
#define SVM_OPT_D2CALGOR_SUPPECACHE_FREEGCACHE          433
#define SVM_OPT_D2CALGOR_SUPPECACHE_FULLGCACHE          434
#define SVM_OPT_D2CALGOR_FULLECACHE_EMPTGCACHE          441
#define SVM_OPT_D2CALGOR_FULLECACHE_DYNAGCACHE          442
#define SVM_OPT_D2CALGOR_FULLECACHE_FREEGCACHE          443
#define SVM_OPT_D2CALGOR_FULLECACHE_FULLGCACHE          444
