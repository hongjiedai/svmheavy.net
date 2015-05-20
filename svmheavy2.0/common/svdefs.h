
/*
 *  SVMheavy - Another SVM Library
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




#ifndef _svdefs_h
#define _svdefs_h

//#define LOOK_SIZE 2061

#include <iostream>
#include <assert.h>
#include "rime.h"

class c_double;

#define INV_ERR         1234321
#define SQRT_ERR        4321234

#define ERR_STREAM std::cerr
#define OUT_STREAM std::cout



//
// This is useful for embedding comments into multiline macros
//

#define MACRO_COMMENT(comment)




//
// ZZERO is the approx range taken to be |zero| (numerical stability)
// ZINF is the approx range taken to be |inf|
//

#define ZZERO 1e-8
#define ZINF 1e30



//
// Function entry recording (if DBGFNENT defined)
//

#ifdef TOP_LEVEL
long _glob_fn_depth_ = 0;
int _glob_fn_rec_on_ = 0;
#endif

#ifndef TOP_LEVEL
extern long _glob_fn_depth_;
extern int _glob_fn_rec_on_;
#endif

#ifndef DBGFNENT
#define FN_ENTRY_POINT
#define FN_EXIT_POINT return
#define FN_RECORD_START
#define FN_RECORD_STOP
#endif

#ifdef DBGFNENT
#define FN_ENTRY_POINT                                                                          \
if ( _glob_fn_rec_on_ )                                                                         \
{                                                                                               \
    long __i___;                                                                                \
    _glob_fn_depth_++;                                                                          \
    for ( __i___ = 1 ; __i___ <= _glob_fn_depth_ ; __i___++ )                                   \
    {                                                                                           \
        std::cerr << " ";                                                                       \
    }                                                                                           \
    std::cerr << "Enter(" << _glob_fn_depth_ << "): " << __FILE__ << "   " << __LINE__ << std::endl; \
}
#define FN_EXIT_POINT                                                                           \
if ( _glob_fn_rec_on_ )                                                                         \
{                                                                                               \
    long __i___;                                                                                \
    for ( __i___ = 1 ; __i___ <= _glob_fn_depth_ ; __i___++ )                                   \
    {                                                                                           \
        std::cerr << " ";                                                                       \
    }                                                                                           \
    std::cerr << "*Exit(" << _glob_fn_depth_ << "): " << __FILE__ << "   " << __LINE__ << std::endl; \
    _glob_fn_depth_--;                                                                          \
}                                                                                               \
return
#define FN_RECORD_START _glob_fn_rec_on_ = 1;
#define FN_RECORD_STOP  _glob_fn_rec_on_ = 0;
#endif




//
// Flop counting terms.
//
// L_DOUBLE - this is used instead of double throughout the code.  If
// flops are not being counted, it defaults to double, so operation is
// as would be expected.  However, if flop counting is done (enabled by
// defining DO__FLOPS), the type c_double is used instead, which acts
// essentially the same as double, except that every floating point
// operation will be counted.
//
// Flop counts are split into kernel and non-kernel using different
// macros defined below.
//
// The counters are as follows:
//
// Flop counters: f_plus  - number of floating point + operations
//                f_minus - number of floating point - operations
//                f_mult  - number of floating point * operations
//                f_div   - number of floating point / operations
//                f_sqrt  - number of floating point sqrt operations
//
// Operations counters: n_iterations - number of iterations
//                      n_step_nonsing - numer of non-singular steps
//                      n_step_sing - numer of singular steps
//                      n_c_lower - variables constrained to lower bound
//                      n_c_upper - variables constrained to upper bound
//                      n_u_lower - variables unconstrained from lower bound
//                      n_u_upper - variables unconstrained from upper bound
//
// Kernel flops: kernel_calls   - number of calls to kernel_dicto fn kernel
//               c_kernel_calls - number of calls to kernel_dicto fn cross_kernel
//               d_kernel_calls - number of calls to kernel_dicto fn dkernel_dkvar
//

#ifndef DO__FLOPS
#define L_DOUBLE double
#endif

#ifdef DO__FLOPS
#define L_DOUBLE c_double
#endif

#ifdef DO__FLOPS

class flop_counts
{
    public:

    unsigned long f_fpu; 
    unsigned long f_plus;
    unsigned long f_minus;
    unsigned long f_mult;
    unsigned long f_div;
    unsigned long f_sqrt;
    unsigned long n_iterations;
    unsigned long n_step_nonsing;
    unsigned long n_step_sing;
    unsigned long n_c_lower;
    unsigned long n_u_lower;
    unsigned long n_c_upper;
    unsigned long n_u_upper;
    unsigned long kernel_calls;
    unsigned long c_kernel_calls;
    unsigned long d_kernel_calls;
};

#define MISC_FLOPS      0
#define SETUP_FLOPS     1
#define OPTIM_FLOPS     2

#ifndef TOP_LEVEL
extern int flop_mode;
extern int do_flop_count;
extern flop_counts flop_counters[3];
#endif

#ifdef TOP_LEVEL
int flop_mode = MISC_FLOPS;
int do_flop_count = 1;
flop_counts flop_counters[3];
#endif

#define STOP_FLOP_COUNT    do_flop_count = 0;
#define RESTART_FLOP_COUNT do_flop_count = 1;

#define SET_FLOP_MODE_MISC      flop_mode = MISC_FLOPS;
#define SET_FLOP_MODE_SETUP     flop_mode = SETUP_FLOPS;
#define SET_FLOP_MODE_OPTIM     flop_mode = OPTIM_FLOPS;

#define REG_f_fpu           if ( do_flop_count ) { ((flop_counters[flop_mode]).f_fpu)++;          }
#define REG_f_plus          if ( do_flop_count ) { ((flop_counters[flop_mode]).f_plus)++;         }
#define REG_f_minus         if ( do_flop_count ) { ((flop_counters[flop_mode]).f_minus)++;        }
#define REG_f_mult          if ( do_flop_count ) { ((flop_counters[flop_mode]).f_mult)++;         }
#define REG_f_div           if ( do_flop_count ) { ((flop_counters[flop_mode]).f_div)++;          }
#define REG_f_sqrt          if ( do_flop_count ) { ((flop_counters[flop_mode]).f_sqrt)++;         }
#define REG_n_iterations    if ( do_flop_count ) { ((flop_counters[flop_mode]).n_iterations)++;   }
#define REG_n_step_nonsing  if ( do_flop_count ) { ((flop_counters[flop_mode]).n_step_nonsing)++; }
#define REG_n_step_sing     if ( do_flop_count ) { ((flop_counters[flop_mode]).n_step_sing)++;    }
#define REG_n_c_lower       if ( do_flop_count ) { ((flop_counters[flop_mode]).n_c_lower)++;      }
#define REG_n_u_lower       if ( do_flop_count ) { ((flop_counters[flop_mode]).n_u_lower)++;      }
#define REG_n_c_upper       if ( do_flop_count ) { ((flop_counters[flop_mode]).n_c_upper)++;      }
#define REG_n_u_upper       if ( do_flop_count ) { ((flop_counters[flop_mode]).n_u_upper)++;      }
#define REG_kernel_calls    if ( do_flop_count ) { ((flop_counters[flop_mode]).kernel_calls)++;   }
#define REG_c_kernel_calls  if ( do_flop_count ) { ((flop_counters[flop_mode]).c_kernel_calls)++; }
#define REG_d_kernel_calls  if ( do_flop_count ) { ((flop_counters[flop_mode]).d_kernel_calls)++; }


#define NON_KERNEL_FLOP_TOTAL (                                 \
(flop_counters[flop_mode]).f_fpu   +                            \
(flop_counters[flop_mode]).f_plus  +                            \
(flop_counters[flop_mode]).f_minus +                            \
(flop_counters[flop_mode]).f_mult  +                            \
(flop_counters[flop_mode]).f_div   +                            \
(flop_counters[flop_mode]).f_sqrt )





#define RESET_FLOPS                                             \
{                                                               \
for ( flop_mode = 0 ; flop_mode <= 2 ; flop_mode++ )            \
{                                                               \
(flop_counters[flop_mode]).f_fpu          = 0;                  \
(flop_counters[flop_mode]).f_plus         = 0;                  \
(flop_counters[flop_mode]).f_minus        = 0;                  \
(flop_counters[flop_mode]).f_mult         = 0;                  \
(flop_counters[flop_mode]).f_div          = 0;                  \
(flop_counters[flop_mode]).f_sqrt         = 0;                  \
(flop_counters[flop_mode]).n_iterations   = 0;                  \
(flop_counters[flop_mode]).n_step_nonsing = 0;                  \
(flop_counters[flop_mode]).n_step_sing    = 0;                  \
(flop_counters[flop_mode]).n_c_lower      = 0;                  \
(flop_counters[flop_mode]).n_u_lower      = 0;                  \
(flop_counters[flop_mode]).n_c_upper      = 0;                  \
(flop_counters[flop_mode]).n_u_upper      = 0;                  \
(flop_counters[flop_mode]).kernel_calls   = 0;                  \
(flop_counters[flop_mode]).c_kernel_calls = 0;                  \
(flop_counters[flop_mode]).d_kernel_calls = 0;                  \
}                                                               \
flop_mode = MISC_FLOPS;                                         \
}


#define PRINT_FLOPS(output)                                                             \
{                                                                                       \
output << "------------------------------------------------\n";                         \
for ( flop_mode = 0 ; flop_mode <= 2 ; flop_mode++ )                                    \
{                                                                                       \
if ( flop_mode == MISC_FLOPS  ) { output << "Miscellaneous flops.\n"; }                 \
if ( flop_mode == SETUP_FLOPS ) { output << "Setup flops.\n";         }                 \
if ( flop_mode == OPTIM_FLOPS ) { output << "Optimisation flops.\n";  }                 \
output << "n_iterations:     " << (flop_counters[flop_mode]).n_iterations   << "\n";    \
output << "n_step_nonsing:   " << (flop_counters[flop_mode]).n_step_nonsing << "\n";    \
output << "n_step_sing:      " << (flop_counters[flop_mode]).n_step_sing    << "\n";    \
output << "n_c_lower:        " << (flop_counters[flop_mode]).n_c_lower      << "\n";    \
output << "n_u_lower:        " << (flop_counters[flop_mode]).n_u_lower      << "\n";    \
output << "n_c_upper:        " << (flop_counters[flop_mode]).n_c_upper      << "\n";    \
output << "n_u_upper:        " << (flop_counters[flop_mode]).n_u_upper      << "\n";    \
output << "kernel_calls:     " << (flop_counters[flop_mode]).kernel_calls   << "\n";    \
output << "c_kernel_calls:   " << (flop_counters[flop_mode]).c_kernel_calls << "\n";    \
output << "d_kernel_calls:   " << (flop_counters[flop_mode]).d_kernel_calls << "\n";    \
output << "Non-Kernel flops: " << NON_KERNEL_FLOP_TOTAL                     << "\n";    \
output << "f_fpu:            " << (flop_counters[flop_mode]).f_fpu          << "\n";    \
output << "f_plus:           " << (flop_counters[flop_mode]).f_plus         << "\n";    \
output << "f_minus:          " << (flop_counters[flop_mode]).f_minus        << "\n";    \
output << "f_mult:           " << (flop_counters[flop_mode]).f_mult         << "\n";    \
output << "f_div:            " << (flop_counters[flop_mode]).f_div          << "\n";    \
output << "f_sqrt:           " << (flop_counters[flop_mode]).f_sqrt         << "\n";    \
}                                                                                       \
flop_mode = MISC_FLOPS;                                                                 \
}


#endif

#ifndef DO__FLOPS

#define STOP_FLOP_COUNT
#define RESTART_FLOP_COUNT
#define NON_KERNEL_FLOP_TOTAL
#define RESET_FLOPS
#define PRINT_FLOPS(output)
#define START_KERNEL_OP
#define END_KERNEL_OP

#endif












//
// Memory debugging - replace occurences of new and delete.
//
// a = new b    should be replaced by REPNEW(a,b);
// a = new b[c] should be replaced by REPNEWB(a,b,c);
//
// delete a   should be replaced be REPDEL(a);
// delete a[] should be replaced be REPDELB(b);
//

#define SAFETY_BUFFER 0
#define PAUSE_LEN -1


#ifdef MEM_DEBUG

#ifndef TOP_LEVEL

extern void **MEM_currents;      // list of current memory assignments
extern char **MEM_currdecs;      // list of descriptors for above
extern int   *MEM_currline;      // list of line numbers for above
extern long MEM_num_currents;    // number of elements in list

extern void **MEM_badrems;       // list of bad removals
extern char **MEM_baddecs;       // list of descriptors for above
extern int   *MEM_badline;       // list of line numbers for above
extern long MEM_num_badrems;     // number of elements in list

extern long MEM_rec_size_curr;   // size of MEM_curr* arrays
extern long MEM_rec_size_badr;   // size of MEM_bad* arrays

extern void **MEM_currents_b;      // list of current memory assignments
extern char **MEM_currdecs_b;      // list of descriptors for above
extern int   *MEM_currline_b;      // list of line numbers for above
extern long MEM_num_currents_b;    // number of elements in list

extern void **MEM_badrems_b;       // list of bad removals
extern char **MEM_baddecs_b;       // list of descriptors for above
extern int   *MEM_badline_b;       // list of line numbers for above
extern long MEM_num_badrems_b;     // number of elements in list

extern long MEM_rec_size_curr_b;   // size of MEM_curr_b* arrays
extern long MEM_rec_size_badr_b;   // size of MEM_bad_b* arrays

#endif

#define MEM_record_ahead_buffer 100

#define MEM_DEFS()                    \
                                      \
        long MEM_rec_size_curr = 0;   \
        long MEM_rec_size_badr = 0;   \
                                      \
        void **MEM_currents = NULL;   \
        char **MEM_currdecs = NULL;   \
        int   *MEM_currline = NULL;   \
        long MEM_num_currents = 0;    \
        void **MEM_badrems = NULL;    \
        char **MEM_baddecs = NULL;    \
        int   *MEM_badline = NULL;    \
        long MEM_num_badrems = 0;     \
                                      \
        long MEM_rec_size_curr_b = 0; \
        long MEM_rec_size_badr_b = 0; \
                                      \
        void **MEM_currents_b = NULL; \
        char **MEM_currdecs_b = NULL; \
        int   *MEM_currline_b = NULL; \
        long MEM_num_currents_b = 0;  \
        void **MEM_badrems_b = NULL;  \
        char **MEM_baddecs_b = NULL;  \
        int   *MEM_badline_b = NULL;  \
        long MEM_num_badrems_b = 0;

#define MEM_REPORT()                                                    \
        {                                                               \
        long ___i,___j;                                                 \
        if ( MEM_num_currents > 0 )                                     \
        {                                                               \
            ___j = 1;                                                   \
            for ( ___i = 1 ; ___i <= MEM_num_currents ; ___i++ )        \
            {                                                           \
                ERR_STREAM << MEM_currline[___i-1] << ": ";             \
                ERR_STREAM << "alloc " << MEM_currdecs[___i-1] << " ";  \
                ERR_STREAM << MEM_currents[___i-1] << "\n";             \
                ___j++;                                                 \
                if ( ___j == PAUSE_LEN )                                \
                {                                                       \
                    std::cin >> ___j;                                   \
                    ___j = 1;                                           \
                }                                                       \
            }                                                           \
        }                                                               \
        if ( MEM_num_badrems > 0 )                                      \
        {                                                               \
            for ( ___i = 1 ; ___i <= MEM_num_badrems ; ___i++ )         \
            {                                                           \
                ERR_STREAM << MEM_badline[___i-1] << ": ";              \
                ERR_STREAM << "delete " << MEM_baddecs[___i-1] << " ";  \
                ERR_STREAM << MEM_badrems[___i-1] << "\n";              \
            }                                                           \
        }                                                               \
        }                                                               \
        {                                                               \
        long ___i,___j;                                                 \
        if ( MEM_num_currents_b > 0 )                                   \
        {                                                               \
            ___j = 1;                                                   \
            for ( ___i = 1 ; ___i <= MEM_num_currents_b ; ___i++ )      \
            {                                                           \
                ERR_STREAM << MEM_currline_b[___i-1] << ": ";           \
                ERR_STREAM << "alloc[] " << MEM_currdecs_b[___i-1] << " "; \
                ERR_STREAM << MEM_currents_b[___i-1] << "\n";           \
                ___j++;                                                 \
                if ( ___j == PAUSE_LEN )                                \
                {                                                       \
                    std::cin >> ___j;                                   \
                    ___j = 1;                                           \
                }                                                       \
            }                                                           \
        }                                                               \
        if ( MEM_num_badrems_b > 0 )                                    \
        {                                                               \
            for ( ___i = 1 ; ___i <= MEM_num_badrems_b ; ___i++ )       \
            {                                                           \
                ERR_STREAM << MEM_badline_b[___i-1] << ": ";            \
                ERR_STREAM << "delete[] " << MEM_baddecs_b[___i-1] << " "; \
                ERR_STREAM << MEM_badrems_b[___i-1] << "\n";            \
            }                                                           \
        }                                                               \
        }

#define REPNEW(name,type)                                               \
        {                                                               \
        name = new type;                                                \
        void **MEM_tempa;                                               \
        char **MEM_tempb;                                               \
        int   *MEM_tempc;                                               \
        long ___i;                                                      \
        if ( ( MEM_rec_size_curr == MEM_num_currents ) &&               \
             ( MEM_rec_size_curr != 0                )    )             \
        {                                                               \
            MEM_rec_size_curr += MEM_record_ahead_buffer;               \
            MEM_tempa = new void *[MEM_rec_size_curr];                  \
            MEM_tempb = new char *[MEM_rec_size_curr];                  \
            MEM_tempc = new int[MEM_rec_size_curr];                     \
            for ( ___i = 1 ; ___i <= MEM_num_currents ; ___i++ )        \
            {                                                           \
                MEM_tempa[___i-1] = MEM_currents[___i-1];               \
                MEM_tempb[___i-1] = MEM_currdecs[___i-1];               \
                MEM_tempc[___i-1] = MEM_currline[___i-1];               \
            }                                                           \
            delete[] MEM_currents;                                      \
            delete[] MEM_currline;                                      \
            delete[] MEM_currdecs;                                      \
            MEM_currents = MEM_tempa;                                   \
            MEM_currdecs = MEM_tempb;                                   \
            MEM_currline = MEM_tempc;                                   \
        }                                                               \
        if ( MEM_rec_size_curr == 0 )                                   \
        {                                                               \
            MEM_rec_size_curr = MEM_record_ahead_buffer;                \
            MEM_currents = new void *[MEM_rec_size_curr];               \
            MEM_currdecs = new char *[MEM_rec_size_curr];               \
            MEM_currline = new int[MEM_rec_size_curr];                  \
        }                                                               \
        MEM_num_currents++;                                             \
        MEM_currents[MEM_num_currents-1] = (void *) name;               \
        MEM_currline[MEM_num_currents-1] = __LINE__;                    \
        MEM_currdecs[MEM_num_currents-1] = __FILE__;                    \
        }

#define REPNEWB(name,type,size)                                         \
        {                                                               \
        name = new type[size+SAFETY_BUFFER];                            \
        void **MEM_tempa;                                               \
        char **MEM_tempb;                                               \
        int   *MEM_tempc;                                               \
        long ___i;                                                      \
        if ( ( MEM_rec_size_curr_b == MEM_num_currents_b ) &&           \
             ( MEM_rec_size_curr_b != 0                  )    )         \
        {                                                               \
            MEM_rec_size_curr_b += MEM_record_ahead_buffer;             \
            MEM_tempa = new void *[MEM_rec_size_curr_b];                \
            MEM_tempb = new char *[MEM_rec_size_curr_b];                \
            MEM_tempc = new int[MEM_rec_size_curr_b];                   \
            for ( ___i = 1 ; ___i <= MEM_num_currents_b ; ___i++ )      \
            {                                                           \
                MEM_tempa[___i-1] = MEM_currents_b[___i-1];             \
                MEM_tempb[___i-1] = MEM_currdecs_b[___i-1];             \
                MEM_tempc[___i-1] = MEM_currline_b[___i-1];             \
            }                                                           \
            delete[] MEM_currents_b;                                    \
            delete[] MEM_currline_b;                                    \
            delete[] MEM_currdecs_b;                                    \
            MEM_currents_b = MEM_tempa;                                 \
            MEM_currdecs_b = MEM_tempb;                                 \
            MEM_currline_b = MEM_tempc;                                 \
        }                                                               \
        if ( MEM_rec_size_curr_b == 0 )                                 \
        {                                                               \
            MEM_rec_size_curr_b = MEM_record_ahead_buffer;              \
            MEM_currents_b = new void *[MEM_rec_size_curr_b];           \
            MEM_currdecs_b = new char *[MEM_rec_size_curr_b];           \
            MEM_currline_b = new int[MEM_rec_size_curr_b];              \
        }                                                               \
        MEM_num_currents_b++;                                           \
        MEM_currents_b[MEM_num_currents_b-1] = (void *) name;           \
        MEM_currline_b[MEM_num_currents_b-1] = __LINE__;                \
        MEM_currdecs_b[MEM_num_currents_b-1] = __FILE__;                \
        }

#define ZREPNEW(name,type)                                              \
        name = new type;

#define ZREPNEWB(name,type,size)                                        \
        name = new type[size+SAFETY_BUFFER];

#define REPDEL(name)                                                    \
        {                                                               \
        long ___i;                                                      \
        long ___j;                                                      \
        int ___was_bad = 0;                                             \
        void *___temp;                                                  \
        ___temp = (void *) name;                                        \
        if ( MEM_num_currents == 0 )                                    \
        {                                                               \
            void **MEM_tempa;                                           \
            char **MEM_tempb;                                           \
            int   *MEM_tempc;                                           \
            if ( ( MEM_rec_size_badr == MEM_num_badrems ) &&            \
                 ( MEM_rec_size_badr != 0               )    )          \
            {                                                           \
                MEM_rec_size_badr += MEM_record_ahead_buffer;           \
                MEM_tempa = new void *[MEM_rec_size_badr];              \
                MEM_tempb = new char *[MEM_rec_size_badr];              \
                MEM_tempc = new int[MEM_rec_size_badr];                 \
                for ( ___i = 1 ; ___i <= MEM_num_badrems ; ___i++ )     \
                {                                                       \
                    MEM_tempa[___i-1] = MEM_badrems[___i-1];            \
                    MEM_tempb[___i-1] = MEM_baddecs[___i-1];            \
                    MEM_tempc[___i-1] = MEM_badline[___i-1];            \
                }                                                       \
                delete[] MEM_badrems;                                   \
                delete[] MEM_badline;                                   \
                delete[] MEM_baddecs;                                   \
                MEM_badrems = MEM_tempa;                                \
                MEM_baddecs = MEM_tempb;                                \
                MEM_badline = MEM_tempc;                                \
            }                                                           \
            if ( MEM_rec_size_badr == 0 )                               \
            {                                                           \
                MEM_rec_size_badr = MEM_record_ahead_buffer;            \
                MEM_badrems = new void *[MEM_rec_size_badr];            \
                MEM_baddecs = new char *[MEM_rec_size_badr];            \
                MEM_badline = new int[MEM_rec_size_badr];               \
            }                                                           \
            MEM_num_badrems++;                                          \
            MEM_badrems[MEM_num_badrems-1] = ___temp;                   \
            MEM_badline[MEM_num_badrems-1] = __LINE__;                  \
            MEM_baddecs[MEM_num_badrems-1] = __FILE__;                  \
            ___was_bad = 1;                                             \
        }                                                               \
        else                                                            \
        {                                                               \
            ___j = 0;                                                   \
            for ( ___i = 1 ; ___i <= MEM_num_currents ; ___i++ )        \
            {                                                           \
                if ( MEM_currents[___i-1] == ___temp )                  \
                {                                                       \
                    ___j = ___i;                                        \
                }                                                       \
            }                                                           \
            if ( ___j == 0 )                                            \
            {                                                           \
                void **MEM_tempa;                                       \
                char **MEM_tempb;                                       \
                int   *MEM_tempc;                                       \
                if ( ( MEM_rec_size_badr == MEM_num_badrems ) &&        \
                     ( MEM_rec_size_badr != 0               )    )      \
                {                                                       \
                    MEM_rec_size_badr += MEM_record_ahead_buffer;       \
                    MEM_tempa = new void *[MEM_rec_size_badr];          \
                    MEM_tempb = new char *[MEM_rec_size_badr];          \
                    MEM_tempc = new int[MEM_rec_size_badr];             \
                    for ( ___i = 1 ; ___i <= MEM_num_badrems ; ___i++ ) \
                    {                                                   \
                        MEM_tempa[___i-1] = MEM_badrems[___i-1];        \
                        MEM_tempb[___i-1] = MEM_baddecs[___i-1];        \
                        MEM_tempc[___i-1] = MEM_badline[___i-1];        \
                    }                                                   \
                    delete[] MEM_badrems;                               \
                    delete[] MEM_badline;                               \
                    delete[] MEM_baddecs;                               \
                    MEM_badrems = MEM_tempa;                            \
                    MEM_baddecs = MEM_tempb;                            \
                    MEM_badline = MEM_tempc;                            \
                }                                                       \
                if ( MEM_rec_size_badr == 0 )                           \
                {                                                       \
                    MEM_rec_size_badr = MEM_record_ahead_buffer;        \
                    MEM_badrems = new void *[MEM_rec_size_badr];        \
                    MEM_baddecs = new char *[MEM_rec_size_badr];        \
                    MEM_badline = new int[MEM_rec_size_badr];           \
                }                                                       \
                MEM_num_badrems++;                                      \
                MEM_badrems[MEM_num_badrems-1] = ___temp;               \
                MEM_badline[MEM_num_badrems-1] = __LINE__;              \
                MEM_baddecs[MEM_num_badrems-1] = __FILE__;              \
                ___was_bad = 1;                                         \
            }                                                           \
            else                                                        \
            {                                                           \
                MEM_num_currents--;                                     \
                if ( MEM_num_currents > 0 )                             \
                {                                                       \
                    for ( ___i = ___j ; ___i <= MEM_num_currents ; ___i++ ) \
                    {                                                   \
                        MEM_currents[___i-1] = MEM_currents[___i];      \
                        MEM_currline[___i-1] = MEM_currline[___i];      \
                        MEM_currdecs[___i-1] = MEM_currdecs[___i];      \
                    }                                                   \
                }                                                       \
            }                                                           \
        }                                                               \
        if ( ___was_bad == 0 )                                          \
        {                                                               \
        delete name;                                                    \
        }                                                               \
        }

#define REPDELB(name)                                                   \
        {                                                               \
        long ___i;                                                      \
        long ___j;                                                      \
        int ___was_bad = 0;                                             \
        void *___temp;                                                  \
        ___temp = (void *) name;                                        \
        if ( MEM_num_currents_b == 0 )                                  \
        {                                                               \
            void **MEM_tempa;                                           \
            char **MEM_tempb;                                           \
            int   *MEM_tempc;                                           \
            if ( ( MEM_rec_size_badr_b == MEM_num_badrems_b ) &&        \
                 ( MEM_rec_size_badr_b != 0                 )    )      \
            {                                                           \
                MEM_rec_size_badr_b += MEM_record_ahead_buffer;         \
                MEM_tempa = new void *[MEM_rec_size_badr_b];            \
                MEM_tempb = new char *[MEM_rec_size_badr_b];            \
                MEM_tempc = new int[MEM_rec_size_badr_b];               \
                for ( ___i = 1 ; ___i <= MEM_num_badrems_b ; ___i++ )   \
                {                                                       \
                    MEM_tempa[___i-1] = MEM_badrems_b[___i-1];          \
                    MEM_tempb[___i-1] = MEM_baddecs_b[___i-1];          \
                    MEM_tempc[___i-1] = MEM_badline_b[___i-1];          \
                }                                                       \
                delete[] MEM_badrems_b;                                 \
                delete[] MEM_badline_b;                                 \
                delete[] MEM_baddecs_b;                                 \
                MEM_badrems_b = MEM_tempa;                              \
                MEM_baddecs_b = MEM_tempb;                              \
                MEM_badline_b = MEM_tempc;                              \
            }                                                           \
            if ( MEM_rec_size_badr_b == 0 )                             \
            {                                                           \
                MEM_rec_size_badr_b = MEM_record_ahead_buffer;          \
                MEM_badrems_b = new void *[MEM_rec_size_badr_b];        \
                MEM_baddecs_b = new char *[MEM_rec_size_badr_b];        \
                MEM_badline_b = new int[MEM_rec_size_badr_b];           \
            }                                                           \
            MEM_num_badrems_b++;                                        \
            MEM_badrems_b[MEM_num_badrems_b-1] = ___temp;               \
            MEM_badline_b[MEM_num_badrems_b-1] = __LINE__;              \
            MEM_baddecs_b[MEM_num_badrems_b-1] = __FILE__;              \
            ___was_bad = 1;                                             \
        }                                                               \
        else                                                            \
        {                                                               \
            ___j = 0;                                                   \
            for ( ___i = 1 ; ___i <= MEM_num_currents_b ; ___i++ )      \
            {                                                           \
                if ( MEM_currents_b[___i-1] == ___temp )                \
                {                                                       \
                    ___j = ___i;                                        \
                }                                                       \
            }                                                           \
            if ( ___j == 0 )                                            \
            {                                                           \
                void **MEM_tempa;                                       \
                char **MEM_tempb;                                       \
                int   *MEM_tempc;                                       \
                if ( ( MEM_rec_size_badr_b == MEM_num_badrems_b ) &&    \
                     ( MEM_rec_size_badr_b != 0                 )    )  \
                {                                                       \
                    MEM_rec_size_badr_b += MEM_record_ahead_buffer;     \
                    MEM_tempa = new void *[MEM_rec_size_badr_b];        \
                    MEM_tempb = new char *[MEM_rec_size_badr_b];        \
                    MEM_tempc = new int[MEM_rec_size_badr_b];           \
                    for ( ___i = 1 ; ___i <= MEM_num_badrems_b ; ___i++ ) \
                    {                                                   \
                        MEM_tempa[___i-1] = MEM_badrems_b[___i-1];      \
                        MEM_tempb[___i-1] = MEM_baddecs_b[___i-1];      \
                        MEM_tempc[___i-1] = MEM_badline_b[___i-1];      \
                    }                                                   \
                    delete[] MEM_badrems_b;                             \
                    delete[] MEM_badline_b;                             \
                    delete[] MEM_baddecs_b;                             \
                    MEM_badrems_b = MEM_tempa;                          \
                    MEM_baddecs_b = MEM_tempb;                          \
                    MEM_badline_b = MEM_tempc;                          \
                }                                                       \
                if ( MEM_rec_size_badr_b == 0 )                         \
                {                                                       \
                    MEM_rec_size_badr_b = MEM_record_ahead_buffer;      \
                    MEM_badrems_b = new void *[MEM_rec_size_badr_b];    \
                    MEM_baddecs_b = new char *[MEM_rec_size_badr_b];    \
                    MEM_badline_b = new int[MEM_rec_size_badr_b];       \
                }                                                       \
                MEM_num_badrems_b++;                                    \
                MEM_badrems_b[MEM_num_badrems_b-1] = ___temp;           \
                MEM_badline_b[MEM_num_badrems_b-1] = __LINE__;          \
                MEM_baddecs_b[MEM_num_badrems_b-1] = __FILE__;          \
                ___was_bad = 1;                                         \
            }                                                           \
            else                                                        \
            {                                                           \
                MEM_num_currents_b--;                                   \
                if ( MEM_num_currents_b > 0 )                           \
                {                                                       \
                    for ( ___i = ___j ; ___i <= MEM_num_currents_b ; ___i++ ) \
                    {                                                   \
                        MEM_currents_b[___i-1] = MEM_currents_b[___i];  \
                        MEM_currline_b[___i-1] = MEM_currline_b[___i];  \
                        MEM_currdecs_b[___i-1] = MEM_currdecs_b[___i];  \
                    }                                                   \
                }                                                       \
            }                                                           \
        }                                                               \
        if ( ___was_bad == 0 )                                          \
        {                                                               \
        delete[] name;                                                  \
        }                                                               \
        }
#endif

#ifndef MEM_DEBUG

#define MEM_DEFS()
#define MEM_REPORT()

#define REPNEW(name,type)                                               \
        name = new type;
#define REPNEWB(name,type,size)                                         \
        name = new type[size+SAFETY_BUFFER];
#define ZREPNEW(name,type)                                              \
        name = new type;
#define ZREPNEWB(name,type,size)                                        \
        name = new type[size+SAFETY_BUFFER];

#define REPDEL(name)                                                    \
        THROW_ASSERT(name != NULL);                                     \
        delete name;                                                    \
        name = NULL;
#define REPDELB(name)                                                   \
        THROW_ASSERT(name != NULL);                                     \
        delete[] name;                                                  \
        name = NULL;

#endif







//
// sanity flag - debug only if flag set
//

#ifndef TOP_LEVEL
extern int sanity_flag;
#endif
#ifdef TOP_LEVEL
int sanity_flag = 1;
#endif









//
// Error handling debugging
//

#ifndef TOP_LEVEL
extern int throw_flag;
extern long throw_val;
extern unsigned long throw_line;
extern char *throw_file;
extern char *throw_date;
extern char *throw_time;
#endif
#ifdef TOP_LEVEL
int throw_flag = 0;
long throw_val = 0;
unsigned long throw_line = 0;
char *throw_file = NULL;
char *throw_date = NULL;
char *throw_time = NULL;
#endif


#define L_THROW(loc_val)                                        \
{                                                               \
throw_flag = 1;                                                 \
throw_val = loc_val;                                            \
throw_line = __LINE__;                                          \
throw_file = __FILE__;                                          \
throw_date = __DATE__;                                          \
throw_time = __TIME__;                                          \
throw loc_val;                                                  \
}

#define REPORT_THROW(output,x)                                  \
{                                                               \
if ( throw_flag )                                               \
{                                                               \
output << "Throw " << throw_val;                                \
output << " on line " << throw_line;                            \
output << " of file " << throw_file << ".\n";                   \
output << "Compilation: " << throw_date;                        \
output << " - " << throw_time << "\n";                          \
}                                                               \
else                                                            \
{                                                               \
output << "Unknown exception " << x << ".\n";                   \
}                                                               \
}



//
// Throwing assertion (so that gdb can catch failure)
//

#ifndef NDEBUG
#define THROW_ASSERT(cond) if ( cond ) {;} else { L_THROW(0); }
#endif
#ifdef NDEBUG
#define THROW_ASSERT(cond)
#endif

#ifndef NDEBUG
#define THROW_ASSERTB(__isf__,cond) if ( cond ) {;} else { L_THROW(__isf__); }
#endif
#ifdef NDEBUG
#define THROW_ASSERTB(__isf___,cond)
#endif



#endif
