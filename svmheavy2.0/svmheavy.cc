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
// Support vector machine - command line interface
//
// Written by: Alistair Shilton
//             Melbourne University
//


#define TOP_LEVEL

#include <windows.h>
#include <iostream>
#include <fstream>
#include <io.h>
#include "nullstream.h"
#include <gsl/err/gsl_errno.h>
#include <time.h>
#include "common/svdefs.h"
#include "common/outfilt.h"
#include "common/sparsevector.h"
#include "svm_pattern.h"
#include "svm_regress.h"
#ifdef USEALLEGRO
#include <allegro.h>
#endif

#ifndef ENOMEM
#define ENOMEM 0
#endif
#ifndef EINVAL
#define EINVAL 0
#endif
#ifndef EPERM
#define EPERM 0
#endif

MEM_DEFS()

unsigned long threadtime;
int *glob_report_interval;

void interact(SVM_pattern &svc);
void interact(SVM_regress &svr);

#ifdef USEALLEGRO
void periodic_report_svc(void *param);
void periodic_report_svr(void *param);
#endif

#ifndef USEALLEGRO
#define END_OF_MAIN()
#endif

void view_help(void);

void view_help(void)
{
    std::cout << "SVMheavy 2.0: an SVM implementation by Alistair Shilton.                    \n"
              << "============                                                                \n"
              << "                                                                            \n"
              << "Usage:     SVMheavy {options} trainfile                                     \n"
              << "trainfile: file containing training vectors.                                \n"
              << "                                                                            \n"
              << "Each line in the training/testing files contains a single sample with the   \n"
              << "following general format (where values in {} are optional):                 \n"
              << "                                                                            \n"
              << "{>,=,<} y/z {tVAL} {TVAL} {eVAL} {EVAL} x                                   \n"
              << "                                                                            \n"
              << "Classification: - {>,=,<} is ignored.                                       \n"
              << "                - y/z gives classification for point (0 for test unknown).  \n"
              << "                - {tVAL} or {TVAL} sets the empirical risk scale.           \n"
              << "                - {eVAL} or {EVAL} sets the distance to surface scale.      \n"
              << "                - x is the training vector.                                 \n"
              << "                                                                            \n"
              << "Regression: - {>,=,<} defines if this is a lower bound, inequality or upper \n"
              << "              bound constraint.                                             \n"
              << "            - y/z gives the target value for this point.                    \n"
              << "            - {tVAL} sets the empirical risk scale for positive errors.     \n"
              << "            - {TVAL} sets the empirical risk scale for negative errors.     \n"
              << "            - {eVAL} sets the epsilon insensitivity for positive errors.    \n"
              << "            - {EVAL} sets the epsilon insensitivity for negative errors.    \n"
              << "            - x is the training vector.                                     \n"
              << "                                                                            \n"
              << "By default, the format for the x vector is:                                 \n"
              << "                                                                            \n"
              << "<feature1>:<value1> <feature2>:<value2> ... <featureN>:<valueN>             \n"
              << "                                                                            \n"
              << "In a testing/training file, other files may be included (a-la #include in C)\n"
              << "using a line like:                                                          \n"
              << "                                                                            \n"
              << "# newfilename                                                               \n"
              << "                                                                            \n"
              << "Pipes: - Non-interactive output sent to standard error.                     \n"
              << "       - Interactive output and help sent to standard out.                  \n"
              << "       - All other output sent direct to files.                             \n"
              << "                                                                            \n"
              << "General options:                                                            \n"
              << "         -?          - this screen.                                         \n"
              << "         -v {0,1,2}  - verbosity level (default 2).                         \n"
              << "                       0 = minimal (write log to trainfile.log).            \n"
              << "                       1 = normal (include alpha, tau in trainfile.log).    \n"
              << "                       2 = maximal (like normal, but also write svm details \n"
              << "                           to trainfile.svm).                               \n"
              << "         -(          - select alternative x format, namely:                 \n"
              << "                       x1 x2 ... xN                                         \n"
              << "         -)          - enter interactive mode after this operation.         \n"
              << "         -C int      - sets report interval, in seconds (default 2).        \n"
              << "Learning options:                                                           \n"
              << "         -z {c,r}    - select between classification 'c' and regression 'r'.\n"
              << "                       (default classification).                            \n"
              << "         -c float    - C: trade-off between training error and margin       \n"
              << "                       (default 1).                                         \n"
              << "         -c+ float   - C+: C scale factor for positive points (default 1).  \n"
              << "         -c- float   - C-: C scale factor for negative points (default 1).  \n"
              << "         -w float    - epsilon width of tube for regression, boundary       \n"
              << "                       distance for classification (default 1 for           \n"
              << "                       classification, 0.001 for regression).               \n"
              << "         -w+ float   - epsilon scaling factor for +ve errors (default 1).   \n"
              << "         -w- float   - epsilon scaling factor for -ve errors (default 1).   \n"
              << "         -b {0,1}    - if 1, use biased hyperplane (i.e. x*w+b0) instead    \n"
              << "                       of unbiased hyperplane (i.e. x*w0) (default 1        \n"
              << "                       (biased)).                                           \n"
              << "         -B float    - if the hyperplane is unbiased, this is used as a     \n"
              << "                       fixed offset for training (i.e. x*w+B) (default 0).  \n"
              << "         -R {0,1}    - use linear linear 0 or quadratic 1 primal empirical  \n"
              << "                       risk (default 0 (linear)).                           \n"
              << "         -T {0,1}    - if 1, use quadratic tube shrinking (regression SVM   \n"
              << "                       only) (default 0 (no tube shrinking)).               \n"
              << "         -! float    - nu value if quadratic tube shrinking used (0.5).     \n"
              << "Preload options:                                                            \n"
              << "         -@ svm_file - preload the SVM from a file.  In this case, the -(,  \n"
              << "                       -z, -b, -B, -R, -T and all optimisation options will \n"
              << "                       be ignored.                                          \n"
              << "Performance estimation options:                                             \n"
              << "         -x {0,1}    - if 1, compute leave-one-out error (0).               \n"
              << "         -$ testfile - if specified, validate svm using testfile.  In this  \n"
              << "                       file, y = 0 indicates unknown classification.        \n"
              << "Kernel options:                                                             \n"
              << "         -t int      - type of kernel function (default 1):                 \n"
              << "                       0 = linear.                                          \n"
              << "                       1 = polynomial (s x*y+c)^d.                          \n"
              << "                       2 = radial basis function exp(-(||x-y||^2)/gamma).   \n"
              << "                       3 = sigmoid tanh(s x*y + c).                         \n"
              << "         -K ker_file - loads the kernel from the file ker_file, overriding  \n"
              << "                       the -t option.                                       \n"
              << "         -d int      - parameter d in polynomial kernel (default 2).        \n"
              << "         -g float    - parameter gamma in rbf kernel (default 1).           \n"
              << "         -s float    - parameter s in sigmoid/poly kernel (default 1).      \n"
              << "         -r float    - parameter c in sigmoid/poly kernel (default 1).      \n"
              << "         -U i x      - for user defined kernels, set uci = x.               \n"
              << "         -V i x      - for user defined kernels, set uvi = x.               \n"
              << "Output options:                                                             \n"
              << "         -L logfile  - string used to derive logfile and svm descript (if   \n"
              << "                       saved).  These files will be this plus an extension  \n"
              << "                       of .log or .svm.                                     \n"
              << "Optimization options:                                                       \n"
              << "         -%          - do not add training points to the SVM.  The trainfile\n"
              << "                       argument is optional in this case.                   \n"
              << "         -*          - do not optimise the SVM.                             \n"
              << "         -e float    - accuracy required of g(x) in solution (defalt 0.001).\n"
              << "         -# int      - terminate optimization after this many iterations,   \n"
              << "                       even if solution not found (default 10000000).       \n"
              << "         -m int      - size of kernel cache in MB (if dynamic, default 40)  \n"
              << "         -W int      - expected training set size (for caching, default 500)\n"
              << "         -A int      - length of buffer used in d2c optimisation method     \n"
              << "                       (default 10).                                        \n"
              << "         -F          - remove non-support vectors after optimisation.       \n"
              << "         -M int      - optimisation method used (default 442).  The format  \n"
              << "                       used is xyz, where x specifies the optimisation      \n"
              << "                       algorithm and yz the caching options.  Specifically: \n"
              << "                       x = 1 for active set method, inverse factorisation.  \n"
              << "                             (NB: the inverse factorisation is numerically  \n"
              << "                                  unstable and hence should not be used).   \n"
              << "                           2 for active set method, cholesky factorisation. \n"
              << "                           3 for Platt's SMO method.                        \n"
              << "                           4 for Daniel's D2C method.                       \n"
              << "                       y = 1 for no gradient caching.                       \n"
              << "                           2 for nonbound support vector gradient caching.  \n"
              << "                           3 for support vector gradient caching.           \n"
              << "                           4 for full gradient caching.                     \n"
              << "                       z = 1 for no kernel caching.                         \n"
              << "                           2 for dynamic kernel caching.                    \n"
              << "                           3 for nonbound support vector kernel caching.    \n"
              << "                           4 for full kernel cache.                         \n"
              << "                       NB: if x = 1 or 2, y must be 3 or 4.                 \n"
              << "                                                                            \n"
              << "Compiled on " << __DATE__ << " at " << __TIME__ << "\n";

    exit(1);

    return;
}

#define FILENAME_SIZE_MAX       128

class SVM_params
{
    public:

    char trainfile[FILENAME_SIZE_MAX];
    char testfile[FILENAME_SIZE_MAX];
    char svm_file[FILENAME_SIZE_MAX];
    char ker_file[FILENAME_SIZE_MAX];
    char log_file[FILENAME_SIZE_MAX];

    int verbosity;
    int validation;
    int vectorform;
    int gointeract;
    int report_interval;

    char svm_type;
    double C;
    double Cplus;
    double Cminus;
    double W;
    double Wplus;
    double Wminus;
    int variable_bias;
    double B;
    int risk_type;
    int tube_type;
    double nu;

    int isC_default;
    int isCplus_default;
    int isCminus_default;
    int isW_default;
    int isWplus_default;
    int isWminus_default;
    int isKer_default;

    int preload_svm;

    int calc_loo;

    int ker_type; /* -1 indicates file used */
    fVECTOR ker_uc;
    fVECTOR ker_uv;
    double ker_xy_scale;

    int no_add;
    int no_opt;
    double sol_tol;
    int opt_method;
    long max_iter;
    int cleanSVs;

    long memsize;
    long min_rowdim;
    long d2cbufflen;
};

void extract_command(int argc, char **argv, SVM_params &p);

void extract_command(int argc, char **argv, SVM_params &p)
{
    int i,j;

    /*
       Set defaults
    */

    (p.trainfile)[0] = '\0';
    (p.testfile)[0]  = '\0';
    (p.svm_file)[0]  = '\0';
    (p.ker_file)[0]  = '\0';
    (p.log_file)[0]  = '\0';

    p.verbosity  = 2;
    p.validation = 0;
    p.vectorform = 1;
    p.gointeract = 0;

    p.svm_type = 'c';
    p.C      = 1;
    p.Cplus  = 1;
    p.Cminus = 1;
    p.W      = 1;
    p.Wplus  = 1;
    p.Wminus = 1;
    p.variable_bias = 1;
    p.B = 0;
    p.risk_type = 0;
    p.tube_type = 0;
    p.nu = 0.5;

    p.isC_default      = 1;
    p.isCplus_default  = 1;
    p.isCminus_default = 1;
    p.isW_default      = 1;
    p.isWplus_default  = 1;
    p.isWminus_default = 1;
    p.isKer_default    = 1;

    p.preload_svm = 0;

    p.calc_loo = 0;

    p.ker_type = 1;
    (p.ker_uc).make_zero(DO_FORCE);
    (p.ker_uc).make_normal(DO_FORCE);
    (p.ker_uc).pad_vector(20);
    p.ker_uc = 0.0;
    (p.ker_uc)[0] = 1.0; /* c for tanh kernel */
    (p.ker_uc)[1] = 1.0; /* c for polynomial kernel */
    (p.ker_uc)[2] = 2.0; /* d */
    (p.ker_uv).make_zero(DO_FORCE);
    (p.ker_uv).make_normal(DO_FORCE);
    (p.ker_uv).pad_vector(20);
    p.ker_uv = 0.0;
    p.ker_xy_scale = 1.0; /* s and gamma */

    p.no_add = 0;
    p.no_opt = 0;
    p.sol_tol = 0.001;
    p.opt_method = 442;
    p.max_iter = 10000000;
    p.cleanSVs = 0;

    p.memsize    = 40;
    p.min_rowdim = 500;
    p.d2cbufflen = 10;

    p.report_interval = 2;

    /*
       parse options
    */

    for ( i = 1 ; i < argc ; i++ )
    {
        if ( argv[i][0] != '-' )
        {
            break;
        }

        switch ( argv[i][1] )
        {
            case '(': { p.vectorform = 0; break; }
            case '*': { p.no_opt     = 1; break; }
            case '%': { p.no_add     = 1; break; }
            case ')': { p.gointeract = 1; break; }
            case 'F': { p.cleanSVs   = 1; break; }

            case 'C': { i++; p.report_interval = atoi(argv[i]); break; }
            case 'v': { i++; p.verbosity       = atoi(argv[i]); break; }
            case 'z': { i++; p.svm_type        = argv[i][0];    break; }
            case 'b': { i++; p.variable_bias   = atoi(argv[i]); break; }
            case 'B': { i++; p.B               = atof(argv[i]); break; }
            case 'R': { i++; p.risk_type       = atoi(argv[i]); break; }
            case '!': { i++; p.nu              = atof(argv[i]); break; }
            case 'x': { i++; p.calc_loo        = atoi(argv[i]); break; }
            case 'e': { i++; p.sol_tol         = atof(argv[i]); break; }
            case 'M': { i++; p.opt_method      = atoi(argv[i]); break; }
            case '#': { i++; p.max_iter        = atoi(argv[i]); break; }
            case 'm': { i++; p.memsize         = atoi(argv[i]); break; }
            case 'W': { i++; p.min_rowdim      = atoi(argv[i]); break; }
            case 'A': { i++; p.d2cbufflen      = atoi(argv[i]); break; }

            case 't': { i++; p.ker_type      = atoi(argv[i]); p.isKer_default = 0; break; }
            case 'd': { i++; (p.ker_uc)[2]   = atof(argv[i]); p.isKer_default = 0; break; }
            case 'g': { i++; p.ker_xy_scale  = atof(argv[i]); p.isKer_default = 0; break; }
            case 's': { i++; p.ker_xy_scale  = atof(argv[i]); p.isKer_default = 0; break; }
            case 'r': { i++; (p.ker_uc)[0]   = atof(argv[i]); p.isKer_default = 0; break; }

            case 'c':
            {
                switch ( argv[i][2] )
                {
                    case '+': { i++; p.Cplus  = atof(argv[i]); p.isCplus_default  = 0; break; }
                    case '-': { i++; p.Cminus = atof(argv[i]); p.isCminus_default = 0; break; }
                    default:  { i++; p.C      = atof(argv[i]); p.isC_default      = 0; break; }
                }

                break;
            }

            case 'w':
            {
                switch ( argv[i][2] )
                {
                    case '+': { i++; p.Wplus  = atof(argv[i]); p.isWplus_default  = 0; break; }
                    case '-': { i++; p.Wminus = atof(argv[i]); p.isWminus_default = 0; break; }
                    default:  { i++; p.W      = atof(argv[i]); p.isW_default      = 0; break; }
                }

                break;
            }

            case 'U':
            {
                p.isKer_default = 0;
                i++;
                j = atoi(argv[i]);
                i++;
                (p.ker_uc)[j-1] = atof(argv[i]);

                break;
            }

            case 'V':
            {
                p.isKer_default = 0;
                i++;
                j = atoi(argv[i]);
                i++;
                (p.ker_uv)[j-1] = atof(argv[i]);
                 
                break;
            }

            case '$':
            {
                p.validation = 1;
                i++;
                strcpy(p.testfile,argv[i]);
                
                break;
            }

            case '@':
            {
                p.preload_svm = 1;
                i++;
                strcpy(p.svm_file,argv[i]);
                
                break;
            }

            case 'L':
            {
                i++;
                strcpy(p.log_file,argv[i]);
                
                break;
            }

            case 'K':
            {
                p.isKer_default = 0;
                p.ker_type = -1;
                i++;
                strcpy(p.ker_file,argv[i]);
                
                break;
            }

            default:
            {
                view_help();

                break;
            }
        }
    }

    if ( p.ker_type == 1 )
    {
        (p.ker_uc)[1] = (p.ker_uc)[0];
        (p.ker_uc)[0] = 1.0;
    }

    /*
       get filename(s)
    */

    if ( i >= argc )
    {
        if ( p.no_add )
        {
            /*
               The trainfile will be taken from the -@ option if available,
               or "unnamed" used otherwise.
            */

            if ( p.preload_svm )
            {
                strcpy(p.trainfile,p.svm_file);
            }

            else
            {
                strcpy(p.trainfile,"unnamed");
            }
        }

        else
        {
            view_help();
        }
    }

    else
    {
        strcpy(p.trainfile,argv[i]);
    }

    /*
       Set W different for regression
    */

    if ( p.isW_default && p.svm_type == 'r' )
    {
        p.W = 0.001;
    }

    return;
}

class SVM_generic
{
    public:

    char type;

    SVM_pattern *svc;
    SVM_regress *svr;
};

#define BUFFER_SIZE 20000

void get_train_info(char *buffer, int &ztype, double &z, double &t, double &tbar, double &eps, double &epsbar, fVECTOR &x, int vectorform);

long process_patternfile(char *buffer, char *filename, SVM_pattern *dest_svm, int vectorform, long id_offset);
long process_regressfile(char *buffer, char *filename, SVM_regress *dest_svm, int vectorform, long id_offset);
long testwith_patternfile(std::ofstream &report_file, char *buffer, char *filename, SVM_pattern *test_svm, int vectorform, long *correct);
long testwith_regressfile(std::ofstream &report_file, char *buffer, char *filename, SVM_regress *test_svm, int vectorform, double *sse);

#ifdef DO__FLOPS
#define PRINT_START_OPT_TIME                                            \
{                                                                       \
    time_t now;                                                         \
    time(&now);                                                         \
    logfile << "operation started  at " << asctime(localtime(&now));    \
    std::cerr << "operation started  at " << asctime(localtime(&now));  \
    RESET_FLOPS                                                         \
}

#define PRINT_END_OPT_TIME                                              \
{                                                                       \
    time_t now;                                                         \
    time(&now);                                                         \
    logfile << "operation finished at " << asctime(localtime(&now));    \
    std::cerr << "operation finished at " << asctime(localtime(&now));  \
    PRINT_FLOPS(logfile)                                                \
}
#endif

#ifndef DO__FLOPS
#define PRINT_START_OPT_TIME                                            \
{                                                                       \
    time_t now;                                                         \
    time(&now);                                                         \
    logfile << "operation started  at " << asctime(localtime(&now));    \
    std::cerr << "operation started  at " << asctime(localtime(&now));  \
}

#define PRINT_END_OPT_TIME                                              \
{                                                                       \
    time_t now;                                                         \
    time(&now);                                                         \
    logfile << "operation finished at " << asctime(localtime(&now));    \
    std::cerr << "operation finished at " << asctime(localtime(&now));  \
}
#endif

#if SVM_HEAVY
int main(int argc, char **argv)
{
    try
    {
        #ifdef USEALLEGRO
        allegro_init();
        install_timer();
        #endif

        gsl_set_error_handler_off();

        SVM_params p;
        SVM_generic machine;
        kernel_dicto K;
        int count;
        char *buffer;

        REPNEWB(buffer,char,BUFFER_SIZE);

        /*
           Get options
        */
   
        extract_command(argc,argv,p);

        /*
           Open logfile.
        */

        if ( (p.log_file)[0] == '\0' )
        {
            strcpy(p.log_file,p.trainfile);
        }

        strcpy(buffer,p.log_file);
        strcat(buffer,".log");

        std::ofstream logfile(buffer);

        if ( !(logfile.is_open()) )
        {
            std::cerr << "Unable to open logfile\n";
            exit(1);
        }

        logfile << "SVM Training details\n";
        logfile << "====================\n\n";

        logfile << "Training file: " << p.trainfile << "\n";
        logfile << "Testing file:  " << p.testfile  << "\n";
        logfile << "Preload file:  " << p.svm_file  << "\n";
        logfile << "Kernel file:   " << p.ker_file  << "\n\n";
        logfile << "Log file:      " << p.log_file  << "\n\n";

        logfile << "Verbosity:       " << p.verbosity       << "\n";
        logfile << "Validate data?:  " << p.validation      << "\n";
        logfile << "Nonsparse data?: " << p.vectorform      << "\n";
        logfile << "Interactive?:    " << p.gointeract      << "\n";
        logfile << "Report interval: " << p.report_interval << "\n\n";

        logfile << "SVM type:       " << p.svm_type      << "\n";
        logfile << "C:              " << p.C             << "\n";
        logfile << "C+:             " << p.Cplus         << "\n";
        logfile << "C-:             " << p.Cminus        << "\n";
        logfile << "W:              " << p.W             << "\n";
        logfile << "W+:             " << p.Wplus         << "\n";
        logfile << "W-:             " << p.Wminus        << "\n";
        logfile << "Variable bias?: " << p.variable_bias << "\n";
        logfile << "Bias (if var):  " << p.B             << "\n";
        logfile << "Risk type:      " << p.risk_type     << "\n";
        logfile << "Shrink tube?:   " << p.tube_type     << "\n";
        logfile << "nu:             " << p.nu            << "\n\n";

        logfile << "C  default?:     " << p.isC_default      << "\n";
        logfile << "C+ default?:     " << p.isCplus_default  << "\n";
        logfile << "C- default?:     " << p.isCminus_default << "\n";
        logfile << "W  default?:     " << p.isW_default      << "\n";
        logfile << "W+ default?:     " << p.isWplus_default  << "\n";
        logfile << "W- default?:     " << p.isWminus_default << "\n";
        logfile << "kernel default?: " << p.isKer_default    << "\n\n";

        logfile << "Preload SVM?:   " << p.preload_svm << "\n";
        logfile << "Calculate LOO?: " << p.calc_loo    << "\n\n";

        logfile << "Kernel type:        " << p.ker_type     << "\n";
        logfile << "Kernel constants:   " << p.ker_uc       << "\n";
        logfile << "Kernel variables:   " << p.ker_uv       << "\n";
        logfile << "XY scale in kernel: " << p.ker_xy_scale << "\n\n";

        logfile << "Don't add vectors?:  " << p.no_add     << "\n";
        logfile << "Don't optimise?:     " << p.no_opt     << "\n";
        logfile << "Solution tolerance:  " << p.sol_tol    << "\n";
        logfile << "Cache size (MB):     " << p.memsize    << "\n";
        logfile << "Expected training:   " << p.min_rowdim << "\n";
        logfile << "d2c buffer length:   " << p.d2cbufflen << "\n";
        logfile << "Optimisation method: " << p.opt_method << "\n";
        logfile << "Maximum iterations:  " << p.max_iter   << "\n";
        logfile << "Cleanup SVMs?:       " << p.cleanSVs   << "\n\n";

        glob_report_interval = &(p.report_interval);

        /*
           Work out machine type.
        */

        if ( p.preload_svm )
        {
            logfile << "Finding type of preload SVM\n";
            logfile << "===========================\n\n";

            std::ifstream svm_source(p.svm_file);

            if ( !(svm_source.is_open()) )
            {
                std::cerr << "Failed to open preload SVM source file.\n";
                exit(1);
            }

            svm_source >> machine.type;
            svm_source.close();
        }

        else
        {
            machine.type = p.svm_type;
        }

        /*
           Construct the kernel function.
        */

        if ( !(p.preload_svm) || !(p.isKer_default) )
        {
            switch ( p.ker_type )
            {
                case -1:
                {
                    logfile << "Loading kernel definition\n";
                    logfile << "=========================\n\n";

                    std::ifstream kern_source(p.ker_file);

                    if ( !(kern_source.is_open()) )
                    {
                        std::cerr << "Failed to open kernel source file.\n";
                        exit(1);
                    }

                    kern_source >> K;
                    kern_source.close();

                    break;
                }

                case 0:
                {
                    logfile << "Constructing linear kernel\n";
                    logfile << "==========================\n\n";

                    fVECTOR zerovect;
                    fMATRIX identmatrix;

                    zerovect.make_zero(DO_FORCE);
                    identmatrix.make_scalar(DO_FORCE);

                    identmatrix[-1] = 1.0;

                    kernel_dicto Ktemp(1.0,1,p.ker_uc,p.ker_uv,1.0,zerovect,identmatrix,1,1,1,0,p.vectorform);

                    K = Ktemp;

                    break;
                }

                case 1:
                {
                    logfile << "Constructing polynomial kernel\n";
                    logfile << "==============================\n\n";

                    fVECTOR zerovect;
                    fMATRIX identmatrix;

                    zerovect.make_zero(DO_FORCE);
                    identmatrix.make_scalar(DO_FORCE);

                    identmatrix[-1] = 1.0;

                    kernel_dicto Ktemp(1.0,6,p.ker_uc,p.ker_uv,sqrt(p.ker_xy_scale),zerovect,identmatrix,1,1,1,0,p.vectorform);

                    K = Ktemp;

                    break;
                }

                case 2:
                {
                    logfile << "Constructing RBF kernel\n";
                    logfile << "=======================\n\n";

                    fVECTOR zerovect;
                    fMATRIX identmatrix;

                    zerovect.make_zero(DO_FORCE);
                    identmatrix.make_scalar(DO_FORCE);

                    identmatrix[-1] = 1.0;

                    kernel_dicto Ktemp(1.0,20,p.ker_uc,p.ker_uv,1.0/sqrt(p.ker_xy_scale),zerovect,identmatrix,1,1,1,0,p.vectorform);

                    K = Ktemp;

                    break;
                }

                case 3:
                {
                    logfile << "Constructing tanh kernel\n";
                    logfile << "========================\n\n";

                    fVECTOR zerovect;
                    fMATRIX identmatrix;

                    zerovect.make_zero(DO_FORCE);
                    identmatrix.make_scalar(DO_FORCE);

                    identmatrix[-1] = 1.0;

                    kernel_dicto Ktemp(1.0,15,p.ker_uc,p.ker_uv,sqrt(p.ker_xy_scale),zerovect,identmatrix,1,1,1,0,p.vectorform);

                    K = Ktemp;

                    break;
                }

                default:
                {
                    view_help();

                    break;
                }
            }
        }

        /*
           Construct the SVM
        */

        machine.svc = NULL;
        machine.svr = NULL;

        if ( 'c' == machine.type )
        {
            if ( p.preload_svm )
            {
                logfile << "Loading SVC from file\n";
                logfile << "=====================\n\n";

                std::cerr << "Loading SVC from file\n";

                std::ifstream svc_source(p.svm_file);

                if ( !(svc_source.is_open()) )
                {
                    std::cerr << "Failed to open SVC source file.\n";
                    exit(1);
                }

                SVM_pattern temp;

                svc_source >> temp;

                REPNEW(machine.svc,SVM_pattern(temp));

                if ( (temp.get_kernel()).is_sparse() )
                {
                    p.vectorform = 1;
                }

                else
                {
                    p.vectorform = 0;
                }

                logfile << "Adjusting SVC parameters\n";
                logfile << "========================\n\n";

                std::cerr << "Adjusting SVC parameters\n";

                if ( !(p.isC_default)      ) { (*(machine.svc)).set_C(p.C);                  }
                if ( !(p.isCplus_default)  ) { (*(machine.svc)).set_C_plus_scale(p.Cplus);   }
                if ( !(p.isCminus_default) ) { (*(machine.svc)).set_C_minus_scale(p.Cminus); }

                if ( !(p.isW_default)      ) { (*(machine.svc)).set_D(p.W);                  }
                if ( !(p.isWplus_default)  ) { (*(machine.svc)).set_D_plus_scale(p.Wplus);   }
                if ( !(p.isWminus_default) ) { (*(machine.svc)).set_D_minus_scale(p.Wminus); }

                if ( !(p.isKer_default) ) { (*(machine.svc)).set_kernel(K); }

                svc_source.close();
            }

            else
            {
                /*
                   Then construct an SVM
                */

                logfile << "Constructing SVC\n";
                logfile << "================\n\n";

                int fix_bias = 1;

                if ( p.variable_bias ) { fix_bias = 0; }

                REPNEW(machine.svc,SVM_pattern(K,p.risk_type,p.C,p.Cplus,p.Cminus,p.W,p.Wplus,p.Wminus,fix_bias,p.B,p.opt_method,p.max_iter,p.sol_tol,p.memsize,p.min_rowdim,p.d2cbufflen));
            }
        }

        else
        {
            if ( p.preload_svm )
            {
                logfile << "Loading SVR from file\n";
                logfile << "=====================\n\n";

                std::cerr << "Loading SVR from file\n";

                std::ifstream svr_source(p.svm_file);

                if ( !(svr_source.is_open()) )
                {
                    std::cerr << "Failed to open SVR source file.\n";
                    exit(1);
                }

                SVM_regress temp;

                svr_source >> temp;

                REPNEW(machine.svr,SVM_regress(temp));

                if ( (temp.get_kernel()).is_sparse() )
                {
                    p.vectorform = 1;
                }

                else
                {
                    p.vectorform = 0;
                }

                logfile << "Adjusting SVR parameters\n";
                logfile << "========================\n\n";

                std::cerr << "Adjusting SVR parameters\n";

                if ( !(p.isC_default)      ) { (*(machine.svr)).set_C(p.C);                  }
                if ( !(p.isCplus_default)  ) { (*(machine.svr)).set_C_plus_scale(p.Cplus);   }
                if ( !(p.isCminus_default) ) { (*(machine.svr)).set_C_minus_scale(p.Cminus); }

                if ( !(p.isW_default)      ) { (*(machine.svr)).set_E(p.W);                  }
                if ( !(p.isWplus_default)  ) { (*(machine.svr)).set_E_plus_scale(p.Wplus);   }
                if ( !(p.isWminus_default) ) { (*(machine.svr)).set_E_minus_scale(p.Wminus); }

                if ( !(p.isKer_default) ) { (*(machine.svc)).set_kernel(K); }

                svr_source.close();
            }

            else
            {
                /*
                   Then construct an SVM
                */

                logfile << "Constructing SVR\n";
                logfile << "================\n\n";

                int fix_bias = 1;

                if ( p.variable_bias ) { fix_bias = 0; }

                REPNEW(machine.svr,SVM_regress(K,p.risk_type,p.tube_type,p.C,p.Cplus,p.Cminus,p.W,p.Wplus,p.Wminus,p.nu,fix_bias,p.B,p.opt_method,p.max_iter,p.sol_tol,p.memsize,p.min_rowdim,p.d2cbufflen));
            }
        }

        /*
           Add new training points
        */

        if ( !(p.no_add) )
        {
            long id_offset;

            if ( 'c' == machine.type )
            {
                id_offset = ((*(machine.svc)).get_max_id())+1;

                logfile << "Adding training points to SVC\n";
                logfile << "=============================\n\n";

                std::cerr << "Adding training points to SVC\n";

                PRINT_START_OPT_TIME

                count = process_patternfile(buffer,p.trainfile,machine.svc,p.vectorform,id_offset);

                PRINT_END_OPT_TIME

                logfile << "Training points added: " << count << "\n\n";
                std::cerr << "Training points added: " << count << "\n";
            }

            else
            {
                id_offset = ((*(machine.svr)).get_max_id())+1;

                logfile << "Adding training points to SVR\n";
                logfile << "=============================\n\n";

                std::cerr << "Adding training points to SVR\n";

                PRINT_START_OPT_TIME

                count = process_regressfile(buffer,p.trainfile,machine.svr,p.vectorform,id_offset);

                PRINT_END_OPT_TIME

                logfile << "Training points added: " << count << "\n\n";
                std::cerr << "Training points added: " << count << "\n";
            }
        }

        /*
           Train the machine
        */

        if ( !(p.no_opt) )
        {
            int threaderr;

            logfile << "Training the SVM\n";
            logfile << "================\n\n";

            std::cerr << "Training the SVM:\n";

            PRINT_START_OPT_TIME

            if ( 'c' == machine.type )
            {
                threadtime = 0;

                if ( ( threaderr = (*(machine.svc)).back_opt_on() ) )
                {
                    logfile   << "Unable to start thread (" << threaderr << ":" << ENOMEM << "," << EINVAL << "," << EPERM << ") - reverting to foreground optimisation\n";
                    std::cerr << "Unable to start thread (" << threaderr << ":" << ENOMEM << "," << EINVAL << "," << EPERM << ") - reverting to foreground optimisation\n";

                    #ifdef USEALLEGRO
                    install_param_int(periodic_report_svc,(void *) (machine.svc),1000*(p.report_interval));
                    #endif

                    count = (*(machine.svc)).optimise();

                    #ifdef USEALLEGRO
                    remove_param_int(periodic_report_svc,(void *) (machine.svc));
                    #endif
                }

                else
                {
                    while ( (*(machine.svc)).get_iter() )
                    {
                        Sleep((p.report_interval));

                        threadtime += (p.report_interval);

                        #ifdef TEST_THREADS
                        (*(machine.svc)).set_C((*(machine.svc)).get_C());
                        #endif

                        std::cerr << threadtime << " sec:  \t";
                        std::cerr << "Iterations: " << (*(machine.svc)).get_iter();
                        std::cerr << "  \tObjective: " << (*(machine.svc)).get_J() << "\n";

                        logfile << threadtime << " sec:  \t";
                        logfile << "Iterations: " << (*(machine.svc)).get_iter();
                        logfile << "  \tObjective: " << (*(machine.svc)).get_J() << "\n";
                    }

                    count = (*(machine.svc)).get_iter_accum();
                }
            }

            else
            {
                threadtime = 0;

                if ( ( threaderr = (*(machine.svr)).back_opt_on() ) )
                {
                    logfile   << "Unable to start thread (" << threaderr << ":" << ENOMEM << "," << EINVAL << "," << EPERM << ") - reverting to foreground optimisation\n";
                    std::cerr << "Unable to start thread (" << threaderr << ":" << ENOMEM << "," << EINVAL << "," << EPERM << ") - reverting to foreground optimisation\n";

                    #ifdef USEALLEGRO
                    install_param_int(periodic_report_svr,(void *) (machine.svr),1000*(p.report_interval));
                    #endif

                    count = (*(machine.svr)).optimise();

                    #ifdef USEALLEGRO
                    remove_param_int(periodic_report_svr,(void *) (machine.svr));
                    #endif
                }

                else
                {
                    while ( (*(machine.svr)).get_iter() )
                    {
                        Sleep((p.report_interval));

                        threadtime += (p.report_interval);

                        std::cerr << threadtime << " sec:  \t";
                        std::cerr << "Iterations: " << (*(machine.svr)).get_iter();
                        std::cerr << "  \tObjective: " << (*(machine.svr)).get_J() << "\n";

                        logfile << threadtime << " sec:  \t";
                        logfile << "Iterations: " << (*(machine.svr)).get_iter();
                        logfile << "  \tObjective: " << (*(machine.svr)).get_J() << "\n";
                    }

                    count = (*(machine.svr)).get_iter_accum();
                }
            }

            PRINT_END_OPT_TIME

            logfile << "Iterations: " << count << "\n\n";
            std::cerr << "Iterations: " << count << "\n";
        }

        /*
           Remove non-support vectors if required
        */

        if ( p.cleanSVs )
        {
            logfile << "Removing non support vectors\n";
            logfile << "============================\n\n";

            std::cerr << "Removing non support vectors\n";

            PRINT_START_OPT_TIME

            if ( 'c' == machine.type )
            {
                (*(machine.svc)).removeNonSupports();
            }

            else
            {
                (*(machine.svr)).removeNonSupports();
            }

            PRINT_END_OPT_TIME
        }

        /*
           Dump information to log file
        */

        {
            std::cerr << "Summary: N    = " << (*(machine.svc)).get_N()    << "\n";
            std::cerr << "         N_C  = " << ((*(machine.svc)).get_N_Z())+((*(machine.svc)).get_N_L())+((*(machine.svc)).get_N_U()) << "\n";
            std::cerr << "         N_S  = " << ((*(machine.svc)).get_N_L())+((*(machine.svc)).get_N_U())+((*(machine.svc)).get_N_F()) << "\n";
            std::cerr << "         N_Z  = " << (*(machine.svc)).get_N_Z()  << "\n";
            std::cerr << "         N_L  = " << (*(machine.svc)).get_N_L()  << "\n";
            std::cerr << "         N_U  = " << (*(machine.svc)).get_N_U()  << "\n";
            std::cerr << "         N_F  = " << (*(machine.svc)).get_N_F()  << "\n";
            std::cerr << "         N_FP = " << (*(machine.svc)).get_N_FP() << "\n";
            std::cerr << "         N_FN = " << (*(machine.svc)).get_N_FN() << "\n";

            logfile << "SVM Summary\n";
            logfile << "===========\n\n";

            logfile << "N    = " << (*(machine.svc)).get_N()    << "\n";
            logfile << "N_C  = " << ((*(machine.svc)).get_N_Z())+((*(machine.svc)).get_N_L())+((*(machine.svc)).get_N_U()) << "\n";
            logfile << "N_S  = " << ((*(machine.svc)).get_N_L())+((*(machine.svc)).get_N_U())+((*(machine.svc)).get_N_F()) << "\n";
            logfile << "N_Z  = " << (*(machine.svc)).get_N_Z()  << "\n";
            logfile << "N_L  = " << (*(machine.svc)).get_N_L()  << "\n";
            logfile << "N_U  = " << (*(machine.svc)).get_N_U()  << "\n";
            logfile << "N_F  = " << (*(machine.svc)).get_N_F()  << "\n";
            logfile << "N_FP = " << (*(machine.svc)).get_N_FP() << "\n";
            logfile << "N_FN = " << (*(machine.svc)).get_N_FN() << "\n\n";

            if ( p.verbosity >= 1 )
            {
                logfile << "\nalpha: \n" << (*(machine.svc)).get_alpha() << "\n";
                logfile << "\ntau:   \n" << (*(machine.svc)).get_tau()   << "\n";
            }
        }

        /*
           Save SVM if applicable.
        */

        if ( p.verbosity >= 2 )
        {
            logfile << "Saving the SVM\n";
            logfile << "==============\n\n";

            if ( 'c' == machine.type )
            {
                char *bufferx;

                REPNEWB(bufferx,char,BUFFER_SIZE);

                strcpy(bufferx,p.log_file);
                strcat(bufferx,".svm");

                std::ofstream fs(bufferx);

                fs << *(machine.svc);

                fs.close();

                REPDELB(bufferx);
            }

            else
            {
                char *bufferx;

                REPNEWB(bufferx,char,BUFFER_SIZE);

                strcpy(bufferx,p.log_file);
                strcat(bufferx,".svm");

                std::ofstream fs(bufferx);

                fs << *(machine.svr);

                fs.close();

                REPDELB(bufferx);
            }
        }

        /*
           Test validity using LOO if applicable
        */

        if ( p.calc_loo )
        {
            logfile << "Calculating the leave-one-out error\n";
            logfile << "===================================\n\n";

            std::cerr << "Calculating the leave-one-out error\n";

            PRINT_START_OPT_TIME

            if ( 'c' == machine.type )
            {
                long looerr;

                looerr = (*(machine.svc)).calc_loo();

                logfile << "Leave-one-out error: " << looerr << "\n\n";
                std::cerr << "Leave-one-out error: " << looerr << "\n";
            }

            PRINT_END_OPT_TIME
        }

        /*
           Test on testing file if applicable
        */

        if ( p.validation )
        {
            if ( 'c' == machine.type )
            {
                long total = 0;
                long correct = 0;

                logfile << "Validating the SVC\n";
                logfile << "==================\n\n";

                std::cerr << "Validating the SVC\n";

                strcpy(buffer,p.testfile);
                strcat(buffer,".res");

                std::ofstream report_file(buffer);

                if ( !(report_file.is_open()) )
                {
                    std::cerr << "Unable to open results file.\n";
                    exit(1);
                }

                PRINT_START_OPT_TIME

                total = testwith_patternfile(report_file,buffer,p.testfile,machine.svc,p.vectorform,&correct);

                PRINT_END_OPT_TIME

                logfile << "Total testing points:                " << total   << "\n";
                logfile << "Correctly classified testing points: " << correct << "\n\n";

                std::cerr << "Total testing points:                " << total   << "\n";
                std::cerr << "Correctly classified testing points: " << correct << "\n";

                report_file.close();
            }

            else
            {
                long total;
                double sse = 0;

                logfile << "Validating the SVR\n";
                logfile << "==================\n\n";

                std::cerr << "Validating the SVR\n";

                strcpy(buffer,p.testfile);
                strcat(buffer,".res");

                std::ofstream report_file(buffer);

                if ( !(report_file.is_open()) )
                {
                    std::cerr << "Unable to open results file.\n";
                    exit(1);
                }

                PRINT_START_OPT_TIME

                total = testwith_regressfile(report_file,buffer,p.testfile,machine.svr,p.vectorform,&sse);

                PRINT_END_OPT_TIME

                logfile << "Total testing points: " << total             << "\n";
                logfile << "Mean squared error:   " << fixnum(sse/total) << "\n\n";

                std::cerr << "Total testing points: " << total             << "\n";
                std::cerr << "Mean squared error:   " << fixnum(sse/total) << "\n";

                report_file.close();
            }
        }

        /*
           Enter interactive mode if required.
        */

        if ( p.gointeract )
        {

            logfile   << "Entering interactive mode\n";
            std::cerr << "Entering interactive mode\n";
            std::cout << "Entering interactive mode\n\n";
            logfile   << "=========================\n\n";
            std::cerr << "=========================\n\n";

            if ( 'c' == machine.type )
            {
                interact(*(machine.svc));
            }

            else
            {
                interact(*(machine.svr));
            }
        }

        /*
           Destroy the SVM
        */

        logfile << "Deleting the SVM\n";
        logfile << "================\n\n";

        if ( NULL != machine.svc )
        {
            REPDEL(machine.svc);
        }

        if ( NULL != machine.svr )
        {
            REPDEL(machine.svr);
        }

        logfile << "Finished!\n";
        logfile << "=========\n\n";

        logfile.close();

        REPDELB(buffer);

        #ifdef USEALLEGRO
        allegro_exit();
        #endif
    }

    catch ( int x )
    {
        REPORT_THROW(std::cerr,x);
    }

    MEM_REPORT();

    return 0;
}
END_OF_MAIN()
#endif

char *jump_stuff(char *buffer);

long process_patternfile(char *buffer, char *filename, SVM_pattern *dest_svm, int vectorform, long id_offset)
{
    double d;
    double t;
    double s;
    int dummy;
    long count = 0;

    /*
       Important to open the file first, as the filename string may
       be overlaid on the buffer.
    */

    std::ifstream infile(filename);

    if ( !(infile.is_open()) )
    {
        std::cerr << "Training file " << filename << " not found\n";
        exit(1);
    }

    while ( !(infile.eof()) )
    {
        infile.getline(buffer,BUFFER_SIZE);

        if ( strlen(buffer) > 0 )
        {
            if ( buffer[0] == '#' )
            {
                char *baffer;

                baffer = buffer + 1;

                while ( ( baffer[0] == ' ' ) || ( baffer[0] == '\t' ) )
                {
                    if ( baffer[0] == '\0' )
                    {
                        std::cerr << "Training file syntax error\n";
                        exit(1);
                    }

                    baffer++;
                }

                count += process_patternfile(buffer,baffer,dest_svm,vectorform,count+id_offset);
            }

            else
            {
                fVECTOR x('r',0);

                get_train_info(buffer,dummy,d,t,t,s,s,x,vectorform);

                dest_svm->add_point(d,x,t,s,count+id_offset);

                count++;
            }
        }
    }

    infile.close();

    return count;
}

long process_regressfile(char *buffer, char *filename, SVM_regress *dest_svm, int vectorform, long id_offset)
{
    double z;
    double t;
    double tbar;
    double eps;
    double epsbar;
    int ztype;
    long count = 0;

    /*
       Important to open the file first, as the filename string may
       be overlaid on the buffer.
    */

    std::ifstream infile(filename);

    if ( !(infile.is_open()) )
    {
        std::cerr << "Training file not found\n";
        exit(1);
    }

    while ( !(infile.eof()) )
    {
        infile.getline(buffer,BUFFER_SIZE);

        if ( strlen(buffer) > 0 )
        {
            if ( buffer[0] == '#' )
            {
                char *baffer;

                baffer = buffer + 1;

                while ( ( baffer[0] == ' ' ) || ( baffer[0] == '\t' ) )
                {
                    if ( baffer[0] == '\0' )
                    {
                        std::cerr << "Training file syntax error\n";
                        exit(1);
                    }

                    baffer++;
                }

                count += process_regressfile(buffer,baffer,dest_svm,vectorform,count+id_offset);
            }

            else
            {
                fVECTOR x('r',0);

                get_train_info(buffer,ztype,z,t,tbar,eps,epsbar,x,vectorform);

                switch ( ztype )
                {
                    case 1:
                    {
                        dest_svm->add_point_lower_bound(z,x,t,eps,count+id_offset);

                        break;
                    }

                    case 2:
                    {
                        dest_svm->add_point_upper_bound(z,x,tbar,epsbar,count+id_offset);

                        break;
                    }

                    default:
                    {
                        dest_svm->add_point(z,x,t,tbar,eps,epsbar,count+id_offset);

                        break;
                    }
                }

                count++;
            }
        }
    }

    infile.close();

    return 0;
}

long testwith_patternfile(std::ofstream &report_file, char *buffer, char *filename, SVM_pattern *test_svm, int vectorform, long *correct)
{
    double d;
    double d_pred;
    double t;
    double s;
    int dummy;
    long count = 0;

    std::ifstream valid_file(filename);

    if ( !(valid_file.is_open()) )
    {
        std::cerr << "Testing file " << filename << " not found.\n";
        exit(1);
    }

    while ( !(valid_file.eof()) )
    {
        valid_file.getline(buffer,BUFFER_SIZE);

        if ( strlen(buffer) > 0 )
        {
            if ( buffer[0] == '#' )
            {
                char *baffer;

                baffer = buffer + 1;

                while ( ( baffer[0] == ' ' ) || ( baffer[0] == '\t' ) )
                {
                    if ( baffer[0] == '\0' )
                    {
                        std::cerr << "Testing file syntax error\n";
                        exit(1);
                    }

                    baffer++;
                }

                count += testwith_patternfile(report_file,buffer,baffer,test_svm,vectorform,correct);
            }

            else
            {
                fVECTOR x('r',0);

                get_train_info(buffer,dummy,d,t,t,s,s,x,vectorform);

                d_pred = test_svm->test_point(x);

                report_file << fixnum(d_pred) << "\n";

                if ( d*d_pred > 0.0 )
                {
                    (*correct)++;
                }

                count++;
            }
        }
    }

    valid_file.close();

    return count;
}

long testwith_regressfile(std::ofstream &report_file, char *buffer, char *filename, SVM_regress *test_svm, int vectorform, double *sse)
{
    double z;
    double z_pred;
    double t;
    double tbar;
    double eps;
    double epsbar;
    int ztype;
    long count = 0;

    std::ifstream valid_file(filename);

    if ( !(valid_file.is_open()) )
    {
        std::cerr << "Testing file " << filename << " not found.\n";
        exit(1);
    }

    while ( !(valid_file.eof()) )
    {
        valid_file.getline(buffer,BUFFER_SIZE);

        if ( strlen(buffer) > 0 )
        {
            if ( buffer[0] == '#' )
            {
                char *baffer;

                baffer = buffer + 1;

                while ( ( baffer[0] == ' ' ) || ( baffer[0] == '\t' ) )
                {
                    if ( baffer[0] == '\0' )
                    {
                        std::cerr << "Testing file syntax error\n";
                        exit(1);
                    }

                    baffer++;
                }

                count += testwith_regressfile(report_file,buffer,baffer,test_svm,vectorform,sse);
            }

            else
            {
                fVECTOR x('r',0);

                get_train_info(buffer,ztype,z,t,tbar,eps,epsbar,x,vectorform);

                z_pred = test_svm->test_point(x);

                report_file << fixnum(z_pred) << "\n";

                switch ( ztype )
                {
                    case 1:
                    {
                        if ( z_pred < z )
                        {
                            (*sse) += ( ( z - z_pred ) * ( z - z_pred ) );
                        }

                        break;
                    }

                    case 2:
                    {
                        if ( z_pred > z )
                        {
                            (*sse) += ( ( z - z_pred ) * ( z - z_pred ) );
                        }

                        break;
                    }

                    default:
                    {
                        (*sse) += ( ( z - z_pred ) * ( z - z_pred ) );

                        break;
                    }
                }

                count++;
            }
        }
    }

    valid_file.close();

    return count;
}

void get_train_info(char *buffer, int &ztype, double &dz, double &t, double &tbar, double &eps, double &epsbar, fVECTOR &x, int vectorform)
{
    if ( strlen(buffer) > 0 )
    {
        /*
           Get type
        */

        ztype = 0;

        switch ( *buffer )
        {
            case '=': { ztype = 0; buffer = jump_stuff(buffer); break; }
            case '>': { ztype = 1; buffer = jump_stuff(buffer); break; }
            case '<': { ztype = 2; buffer = jump_stuff(buffer); break; }

            default:
            {
                break;
            }
        }

        /*
           First get z
        */
                    
        sscanf(buffer,"%lf",&dz);
        buffer = jump_stuff(buffer);

        /*
           Get additional information (if present)
        */

        t      = 1;
        tbar   = 1;
        eps    = 1;
        epsbar = 1;

        while ( ( *buffer == 't' ) || ( *buffer == 'T' ) ||
                ( *buffer == 'e' ) || ( *buffer == 'E' )    )
        {
            switch ( *buffer )
            {
                case 't': { sscanf(buffer+1,"%lf",&t);      break; }
                case 'T': { sscanf(buffer+1,"%lf",&tbar);   break; }
                case 'e': { sscanf(buffer+1,"%lf",&eps);    break; }
                case 'E': { sscanf(buffer+1,"%lf",&epsbar); break; }

                default:
                {
                    L_THROW(0);

                    break;
                }
            }

            buffer = jump_stuff(buffer); 
        }

        /*
           Get vector and add point.
        */

        if ( vectorform )
        {
            convert_standard_format(buffer,x);
        }

        else
        {
            convert_full_format(buffer,x);
        }
    }

    return;
}

char *jump_stuff(char *buffer)
{
    while ( !isspace(*buffer) )
    {
        buffer++;
    }

    while ( isspace(*buffer) )
    {
        buffer++;
    }

    return buffer;
}

#define XPRINT_START_ADD_TIME                                            \
{                                                                       \
    time_t now;                                                         \
    time(&now);                                                         \
    std::cout << "priming was started at " << asctime(localtime(&now)) << "\n"; \
}

#define XPRINT_END_ADD_TIME                                              \
{                                                                       \
    time_t now;                                                         \
    time(&now);                                                         \
    std::cout << "priming was completed at " << asctime(localtime(&now)) << "\n"; \
}

#ifdef DO__FLOPS
#define XPRINT_START_OPT_TIME                                            \
{                                                                       \
    time_t now;                                                         \
    time(&now);                                                         \
    std::cout << "training was started at " << asctime(localtime(&now)) << "\n"; \
    RESET_FLOPS                                                         \
}

#define XPRINT_END_OPT_TIME                                              \
{                                                                       \
    time_t now;                                                         \
    time(&now);                                                         \
    std::cout << "training was completed at " << asctime(localtime(&now)) << "\n"; \
    PRINT_FLOPS(std::cout)                                              \
}
#endif

#ifndef DO__FLOPS
#define XPRINT_START_OPT_TIME                                            \
{                                                                       \
    time_t now;                                                         \
    time(&now);                                                         \
    std::cout << "training was started at " << asctime(localtime(&now)) << "\n"; \
}

#define XPRINT_END_OPT_TIME                                              \
{                                                                       \
    time_t now;                                                         \
    time(&now);                                                         \
    std::cout << "training was completed at " << asctime(localtime(&now)) << "\n"; \
}
#endif

void interact(SVM_pattern &machine)
{
    {
        int echo = 0;
        long many;
        long i,j;
        char filename[256];
        int id_count = 1;

        std::cout << "Echo (0 = off, 1 = on 2 = quiet echo): ";
        std::cin >> echo;
        
        repeat:

        std::cout << "Options: 1 = Add training point.\n";
        std::cout << "         2 = Add training points.\n";
        std::cout << "         3 = Add training points from file.\n";
        std::cout << "         4 = Add training points from file (no optimise).\n";
        std::cout << "         5 = Add training points from file (sparse).\n";
        std::cout << "         6 = Add training points from file (sparse, no opt).\n";
        std::cout << "         7 = Set C/N.\n";
        std::cout << "         8 = Set C+ scale.\n";
        std::cout << "         9 = Set C- scale.\n";
        std::cout << "         10 = Set D.\n";
        std::cout << "         11 = Set D+ scale.\n";
        std::cout << "         12 = Set D- scale.\n";

        std::cout << "         14 = Set kernel.\n";
        std::cout << "         15 = Delete point.\n";
        std::cout << "         16 = Test points.\n";
        std::cout << "         17 = Test performance.\n";
        std::cout << "         18 = Test points from file.\n";
        std::cout << "         19 = Test performance from file.\n";
        std::cout << "         20 = Test points from file (sparse).\n";
        std::cout << "         21 = Test performance from file (sparse).\n";
        std::cout << "         22 = Print SVM.\n";
        std::cout << "         23 = Save SVM.\n";
        std::cout << "         24 = Load SVM.\n";
        std::cout << "         25 = Veiw report.\n";
        std::cout << "         32 = Optimise SVM.\n";
        std::cout << "         33 = Optimise SVM to specified accuracy.\n";
        std::cout << "         34 = Re-optimise SVM to specified accuracy with time limit.\n";
        std::cout << "         35 = Calculate leave-one-out error.\n";

        std::cin >> many; if ( echo ) { std::cout << many << "\n"; }

        switch ( many )
        {
            case 1:
            {
                L_DOUBLE d;
                fVECTOR x;
                L_DOUBLE t;

                L_DOUBLE s;


                std::cout << "d" << id_count << ": "; std::cin >> d; if ( echo ) { std::cout << d << "\n"; }

                if ( echo == 1 ) { x.prime_io(&(std::cout),1); }
                else             { x.prime_io(&(std::cout),0); }

                std::cout << "x:\n\n"; std::cin >> x;

                std::cout << "t: "; std::cin >> t; if ( echo ) { std::cout << t << "\n"; }

                std::cout << "s: "; std::cin >> s; if ( echo ) { std::cout << s << "\n"; }


                std::cout << "\n";

                machine.add_point(d,x,t,s,id_count);
                machine.optimise();

                id_count++;

                break;
            }

            case 2:
            {
                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }

                fVECTOR d('r',many);
                f_fVECTOR x('r',many);
                fVECTOR t('r',many);

                fVECTOR s('r',many);

                iVECTOR id('r',many);

                for ( i = 1 ; i <= many ; i++ )
                {
                    std::cout << "d" << id_count << ": "; std::cin >> d[i-1]; if ( echo ) { std::cout << d[i-1] << "\n"; }

                    if ( echo == 1 ) { (x[i-1]).prime_io(&(std::cout),1); }
                    else             { (x[i-1]).prime_io(&(std::cout),0); }

                    std::cout << "x:\n\n"; std::cin >> x[i-1];

                    std::cout << "t: "; std::cin >> t[i-1]; if ( echo ) { std::cout << t[i-1] << "\n"; }

                    std::cout << "s: "; std::cin >> s[i-1]; if ( echo ) { std::cout << s[i-1] << "\n"; }


                    std::cout << "\n";

                    id[i-1] = id_count;

                    id_count++;
                }

                XPRINT_START_ADD_TIME
                machine.add_points(d,x,t,s,id);
                XPRINT_END_ADD_TIME

                XPRINT_START_OPT_TIME
                FN_RECORD_START
                machine.optimise();
                FN_RECORD_STOP
                XPRINT_END_OPT_TIME

                break;
            }

            case 3:
            {
                long dim;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "Dimension: "; std::cin >> dim; if ( echo ) { std::cout << dim << "\n"; }
                std::cout << "File name (d x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ifstream gs(filename);

                fVECTOR d('r',many);
                f_fVECTOR x('r',many);
                fVECTOR t('r',many);

                fVECTOR s('r',many);

                iVECTOR id('r',many);

                for ( i = 1 ; i <= many ; i++ )
                {
                    gs >> d[i-1];

                    (x[i-1]).make_normal(DO_FORCE);
                    (x[i-1]).pad_vector(dim);

                    for ( j = 1 ; j <= dim ; j++ )
                    {
                        gs >> (x[i-1])[j-1];
                    }

                    t[i-1] = 1.0;

                    s[i-1] = 1.0;


                    id[i-1] = id_count;

                    id_count++;
                }

                gs.close();

                XPRINT_START_ADD_TIME
                machine.add_points(d,x,t,s,id);
                XPRINT_END_ADD_TIME

                XPRINT_START_OPT_TIME
                FN_RECORD_START
                machine.optimise();
                FN_RECORD_STOP
                XPRINT_END_OPT_TIME

                break;
            }

            case 4:
            {
                long dim;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "Dimension: "; std::cin >> dim; if ( echo ) { std::cout << dim << "\n"; }
                std::cout << "File name (d x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ifstream gs(filename);

                fVECTOR d('r',many);
                f_fVECTOR x('r',many);
                fVECTOR t('r',many);

                fVECTOR s('r',many);

                iVECTOR id('r',many);

                for ( i = 1 ; i <= many ; i++ )
                {
                    gs >> d[i-1];

                    (x[i-1]).make_normal(DO_FORCE);
                    (x[i-1]).pad_vector(dim);

                    for ( j = 1 ; j <= dim ; j++ )
                    {
                        gs >> (x[i-1])[j-1];
                    }

                    t[i-1] = 1.0;

                    s[i-1] = 1.0;


                    id[i-1] = id_count;

                    id_count++;
                }

                gs.close();

                machine.add_points(d,x,t,s,id);

                break;
            }

            case 5:
            {
                long dim;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "File name (d x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ifstream gs(filename);

                fVECTOR d('r',many);
                f_fVECTOR x('r',many);
                fVECTOR t('r',many);

                fVECTOR s('r',many);

                iVECTOR id('r',many);

                for ( i = 1 ; i <= many ; i++ )
                {
                    gs >> d[i-1];

                    gs >> dim;

                    (x[i-1]).make_normal(DO_FORCE);
                    (x[i-1]).pad_vector(dim);

                    for ( j = 1 ; j <= dim ; j++ )
                    {
                        gs >> (x[i-1])[j-1];
                    }

                    t[i-1] = 1.0;

                    s[i-1] = 1.0;


                    id[i-1] = id_count;

                    id_count++;
                }

                gs.close();

                XPRINT_START_ADD_TIME
                machine.add_points(d,x,t,s,id);
                XPRINT_END_ADD_TIME

                XPRINT_START_OPT_TIME
                FN_RECORD_START
                machine.optimise();
                FN_RECORD_STOP
                XPRINT_END_OPT_TIME

                break;
            }

            case 6:
            {
                long dim;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "File name (d x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ifstream gs(filename);

                fVECTOR d('r',many);
                f_fVECTOR x('r',many);
                fVECTOR t('r',many);

                fVECTOR s('r',many);

                iVECTOR id('r',many);

                for ( i = 1 ; i <= many ; i++ )
                {
                    gs >> d[i-1];

                    gs >> dim;

                    (x[i-1]).make_normal(DO_FORCE);
                    (x[i-1]).pad_vector(dim);

                    for ( j = 1 ; j <= dim ; j++ )
                    {
                        gs >> (x[i-1])[j-1];
                    }

                    t[i-1] = 1.0;

                    s[i-1] = 1.0;


                    id[i-1] = id_count;

                    id_count++;
                }

                gs.close();

                machine.add_points(d,x,t,s,id);

                break;
            }

            case 7:
            {
                L_DOUBLE C_new;

                std::cout << "C: "; std::cin >> C_new; if ( echo ) { std::cout << C_new << "\n"; }

                machine.set_C(C_new);
                machine.optimise();

                break;
            }

            case 8:
            {
                L_DOUBLE C_new;

                std::cout << "C+ scale: "; std::cin >> C_new; if ( echo ) { std::cout << C_new << "\n"; }

                machine.set_C_plus_scale(C_new);
                machine.optimise();

                break;
            }

            case 9:
            {
                L_DOUBLE C_new;

                std::cout << "C- scale: "; std::cin >> C_new; if ( echo ) { std::cout << C_new << "\n"; }

                machine.set_C_minus_scale(C_new);
                machine.optimise();

                break;
            }

            case 10:
            {
                L_DOUBLE D_new;

                std::cout << "D: "; std::cin >> D_new; if ( echo ) { std::cout << D_new << "\n"; }

                machine.set_D(D_new);
                machine.optimise();

                break;
            }

            case 11:
            {
                L_DOUBLE D_new;

                std::cout << "D+ scale: "; std::cin >> D_new; if ( echo ) { std::cout << D_new << "\n"; }

                machine.set_D_plus_scale(D_new);
                machine.optimise();

                break;
            }

            case 12:
            {
                L_DOUBLE D_new;

                std::cout << "D- scale: "; std::cin >> D_new; if ( echo ) { std::cout << D_new << "\n"; }

                machine.set_D_minus_scale(D_new);
                machine.optimise();

                break;
            }

            case 14:
            {
                kernel_dicto K;

                std::cout << "Enter kernel details (sparsity must match training data).\n\n";

                K.prime_io(&(std::cout),echo);

                std::cin >> K;

                if ( (machine.get_kernel()).is_sparse() != K.is_sparse() )
                {
                    std::cout << "Sparsity mismatch - kernel unchanged.\n";
                }

                else
                {
                    machine.set_kernel(K);
                    machine.optimise();
                }

                break;
            }

            case 15:
            {
                std::cout << "Which point: "; std::cin >> i; if ( echo ) { std::cout << i << "\n"; }

                machine.del_point(i);
                machine.optimise();

                break;
            }

            case 16:
            {
                fVECTOR x;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }

                for ( i = 1 ; i <= many ; i++ )
                {
                    if ( echo == 1 ) { x.prime_io(&(std::cout),1); }
                    else             { x.prime_io(&(std::cout),0); }

                    std::cout << "x:\n\n"; std::cin >> x;

                    std::cout << "g(x) = " << machine.test_point(x) << "\n\n";
                }

                break;
            }

            case 17:
            {
                fVECTOR x;
                L_DOUBLE d;
                L_DOUBLE g;
                long errors;



                errors = 0;


                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }

                for ( i = 1 ; i <= many ; i++ )
                {
                    std::cout << "d: "; std::cin >> d; if ( echo ) { std::cout << d << "\n"; }

                    if ( echo == 1 ) { x.prime_io(&(std::cout),1); }
                    else             { x.prime_io(&(std::cout),0); }

                    std::cout << "x:\n\n"; std::cin >> x;

                    g = machine.test_point(x);

                    std::cout << "g(x) = " << g << "\n\n";

                    if ( d*g <= 0.0 ) { errors++; }







                }



                std::cout << "Errors: " << errors << "\n\n";



                break;
            }

            case 18:
            {
                long dim;
                long startpoint;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "Number of points to ignore at start: "; std::cin >> startpoint; if ( echo ) { std::cout << startpoint << "\n"; }
                std::cout << "Dimension: "; std::cin >> dim; if ( echo ) { std::cout << dim << "\n"; }
                std::cout << "Input file name (x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ifstream gs(filename);

                std::cout << "Report file name (g(x) format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ofstream hs(filename);

                fVECTOR x('r',dim);

                for ( i = 1 ; i <= many ; i++ )
                {
                    if ( i > startpoint )
                    {
                        for ( j = 1 ; j <= dim ; j++ )
                        {
                            gs >> x[j-1];
                        }
                    }

                    hs << machine.test_point(x) << "\n";
                }

                gs.close();
                hs.close();

                break;
            }

            case 19:
            {
                long dim;
                long startpoint;
                long numpoints;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "Number of points to ignore at start: "; std::cin >> startpoint; if ( echo ) { std::cout << startpoint << "\n"; }
                std::cout << "Dimension: "; std::cin >> dim; if ( echo ) { std::cout << dim << "\n"; }
                std::cout << "Input file name (d x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                numpoints = many-startpoint;

                std::ifstream gs(filename);

                std::cout << "Report file name (g(x) format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ofstream hs(filename);

                fVECTOR x('r',dim);
                long errors;


                L_DOUBLE d;
                L_DOUBLE g;

                errors = 0;


                for ( i = 1 ; i <= many ; i++ )
                {
                    gs >> d;

                    for ( j = 1 ; j <= dim ; j++ )
                    {
                        gs >> x[j-1];
                    }

                    if ( i > startpoint )
                    {
                        g = machine.test_point(x);

                        hs << g << "\n";

                        if ( g*d <= 0.0 ) { errors++; }







                    }
                }



                std::cout << "Errors: " << errors << "\n\n";



                gs.close();
                hs.close();

                break;
            }

            case 20:
            {
                long dim;
                long startpoint;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "Number of points to ignore at start: "; std::cin >> startpoint; if ( echo ) { std::cout << startpoint << "\n"; }
                std::cout << "Input file name (x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ifstream gs(filename);

                std::cout << "Report file name (g(x) format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ofstream hs(filename);

                fVECTOR x('r',dim);

                for ( i = 1 ; i <= many ; i++ )
                {
                    gs >> dim;

                    if ( i > startpoint )
                    {
                        for ( j = 1 ; j <= dim ; j++ )
                        {
                            gs >> x[j-1];
                        }
                    }

                    hs << machine.test_point(x) << "\n";
                }

                gs.close();
                hs.close();

                break;
            }

            case 21:
            {
                long dim;
                long startpoint;
                long numpoints;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "Number of points to ignore at start: "; std::cin >> startpoint; if ( echo ) { std::cout << startpoint << "\n"; }
                std::cout << "Input file name (d x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                numpoints = many-startpoint;

                std::ifstream gs(filename);

                std::cout << "Report file name (g(x) format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ofstream hs(filename);

                fVECTOR x('r',10);
                long errors;


                L_DOUBLE d;
                L_DOUBLE g;

                errors = 0;


                for ( i = 1 ; i <= many ; i++ )
                {
                    gs >> d;

                    gs >> dim;

                    x.make_zero(DO_FORCE);
                    x.make_normal(DO_FORCE);
                    x.pad_vector(dim);

                    for ( j = 1 ; j <= dim ; j++ )
                    {
                        gs >> x[j-1];
                    }

                    if ( i > startpoint )
                    {
                        g = machine.test_point(x);

                        hs << g << "\n";

                        if ( g*d <= 0.0 ) { errors++; }







                    }
                }



                std::cout << "Errors: " << errors << "\n\n";



                gs.close();
                hs.close();

                break;
            }

            case 22:
            {
                std::cout << machine << "\n";

                break;
            }

            case 23:
            {
                std::cout << "File name: ";
                std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ofstream fs(filename);

                fs << machine;

                fs.close();

                break;
            }

            case 24:
            {
                std::cout << "File name: ";
                std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ifstream gs(filename);

                gs >> machine;

                gs.close();

                break;
            }

            case 25:
            {
                std::cout << "Size:                          " << machine.get_N()   << "\n";
                std::cout << "Support vectors:               " << machine.get_N_S() << "\n";
                std::cout << "Unconstrained support vectors: " << machine.get_N_F() << "\n";

                break;
            }

            case 32:
            {
                XPRINT_START_OPT_TIME
                FN_RECORD_START
                machine.optimise();
                FN_RECORD_STOP
                XPRINT_END_OPT_TIME

                break;
            }

            case 33:
            {
                L_DOUBLE local_tol;

                std::cout << "Tolerance: "; std::cin >> local_tol; if ( echo ) { std::cout << local_tol << "\n"; }

                XPRINT_START_OPT_TIME
                FN_RECORD_START
                machine.optimise(local_tol);
                FN_RECORD_STOP
                XPRINT_END_OPT_TIME

                break;
            }

            case 34:
            {
                L_DOUBLE local_tol;
                int max_iter;

                std::cout << "Tolerance: ";      std::cin >> local_tol; if ( echo ) { std::cout << local_tol << "\n"; }
                std::cout << "Max iterations: "; std::cin >> max_iter;  if ( echo ) { std::cout << max_iter  << "\n"; }

                XPRINT_START_OPT_TIME
                FN_RECORD_START
                machine.optimise_from_zero(local_tol,max_iter);
                FN_RECORD_STOP
                XPRINT_END_OPT_TIME

                break;
            }

            case 35:
            {
                std::cout << "Leave-one-out error: " << machine.calc_loo() << "/" << machine.get_N() << "\n";

                break;
            }

            default:
            {
                return;

                break;
            }
        }

        goto repeat;
    }

    return;
}

void interact(SVM_regress &machine)
{
    {
        int echo = 0;
        long many;
        long i,j;
        char filename[256];
        int id_count = 1;

        std::cout << "Echo (0 = off, 1 = on 2 = quiet echo): ";
        std::cin >> echo;

        repeat:

        std::cout << "Options: 1 = Add training point.\n";
        std::cout << "         2 = Add training points.\n";
        std::cout << "         3 = Add training points from file.\n";
        std::cout << "         4 = Add training points from file (no optimise).\n";
        std::cout << "         5 = Add training points from file (sparse).\n";
        std::cout << "         6 = Add training points from file (sparse, no opt).\n";
        std::cout << "         7 = Set C/N.\n";
        std::cout << "         8 = Set C+ scale.\n";
        std::cout << "         9 = Set C- scale.\n";
        std::cout << "         10 = Set E.\n";
        std::cout << "         11 = Set E+ scale.\n";
        std::cout << "         12 = Set E- scale.\n";
        std::cout << "         13 = Set nu.\n";
        std::cout << "         14 = Set kernel.\n";
        std::cout << "         15 = Delete point.\n";
        std::cout << "         16 = Test points.\n";
        std::cout << "         17 = Test performance.\n";
        std::cout << "         18 = Test points from file.\n";
        std::cout << "         19 = Test performance from file.\n";
        std::cout << "         20 = Test points from file (sparse).\n";
        std::cout << "         21 = Test performance from file (sparse).\n";
        std::cout << "         22 = Print SVM.\n";
        std::cout << "         23 = Save SVM.\n";
        std::cout << "         24 = Load SVM.\n";
        std::cout << "         25 = Veiw report.\n";






        std::cout << "         33 = Optimise SVM.\n";
        std::cout << "         34 = Optimise SVM to specified accuracy.\n";
        std::cout << "         35 = Re-optimise SVM to specified accuracy with time limit.\n";


        std::cin >> many; if ( echo ) { std::cout << many << "\n"; }

        switch ( many )
        {
            case 1:
            {
                L_DOUBLE z;
                fVECTOR x;
                L_DOUBLE t;
                L_DOUBLE tbar;
                L_DOUBLE eps;
                L_DOUBLE eps_bar;

                std::cout << "z" << id_count << ": "; std::cin >> z; if ( echo ) { std::cout << z << "\n"; }

                if ( echo == 1 ) { x.prime_io(&(std::cout),1); }
                else             { x.prime_io(&(std::cout),0); }

                std::cout << "x:\n\n"; std::cin >> x;

                std::cout << "t: "; std::cin >> t; if ( echo ) { std::cout << t << "\n"; }
                std::cout << "tbar: "; std::cin >> tbar; if ( echo ) { std::cout << tbar << "\n"; }
                std::cout << "epsilon scale: "; std::cin >> eps; if ( echo ) { std::cout << eps << "\n"; }
                std::cout << "epsilon_bar scale: "; std::cin >> eps_bar; if ( echo ) { std::cout << eps << "\n"; }

                std::cout << "\n";

                machine.add_point(z,x,t,tbar,eps,eps_bar,id_count);
                machine.optimise();

                id_count++;

                break;
            }

            case 2:
            {
                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }

                fVECTOR z('r',many);
                f_fVECTOR x('r',many);
                fVECTOR t('r',many);
                fVECTOR tbar('r',many);
                fVECTOR eps('r',many);
                fVECTOR eps_bar('r',many);
                iVECTOR id('r',many);

                for ( i = 1 ; i <= many ; i++ )
                {
                    std::cout << "z" << id_count << ": "; std::cin >> z[i-1]; if ( echo ) { std::cout << z[i-1] << "\n"; }

                    if ( echo == 1 ) { (x[i-1]).prime_io(&(std::cout),1); }
                    else             { (x[i-1]).prime_io(&(std::cout),0); }

                    std::cout << "x:\n\n"; std::cin >> x[i-1];

                    std::cout << "t: "; std::cin >> t[i-1]; if ( echo ) { std::cout << t[i-1] << "\n"; }
                    std::cout << "tbar: "; std::cin >> tbar[i-1]; if ( echo ) { std::cout << tbar[i-1] << "\n"; }
                    std::cout << "epsilon scale: "; std::cin >> eps[i-1]; if ( echo ) { std::cout << eps[i-1] << "\n"; }
                    std::cout << "epsilon_bar scale: "; std::cin >> eps_bar[i-1]; if ( echo ) { std::cout << eps[i-1] << "\n"; }

                    std::cout << "\n";

                    id[i-1] = id_count;

                    id_count++;
                }

                XPRINT_START_ADD_TIME
                machine.add_points(z,x,t,tbar,eps,eps_bar,id);
                XPRINT_END_ADD_TIME

                XPRINT_START_OPT_TIME
                FN_RECORD_START
                machine.optimise();
                FN_RECORD_STOP
                XPRINT_END_OPT_TIME

                break;
            }

            case 3:
            {
                long dim;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "Dimension: "; std::cin >> dim; if ( echo ) { std::cout << dim << "\n"; }
                std::cout << "File name (z x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ifstream gs(filename);

                fVECTOR z('r',many);
                f_fVECTOR x('r',many);
                fVECTOR t('r',many);
                fVECTOR tbar('r',many);
                fVECTOR eps('r',many);
                fVECTOR eps_bar('r',many);
                iVECTOR id('r',many);

                for ( i = 1 ; i <= many ; i++ )
                {
                    gs >> z[i-1];

                    (x[i-1]).make_normal(DO_FORCE);
                    (x[i-1]).pad_vector(dim);

                    for ( j = 1 ; j <= dim ; j++ )
                    {
                        gs >> (x[i-1])[j-1];
                    }

                    t[i-1] = 1.0;
                    tbar[i-1] = 1.0;
                    eps[i-1] = 1.0;
                    eps_bar[i-1] = 1.0;

                    id[i-1] = id_count;

                    id_count++;
                }

                gs.close();

                XPRINT_START_ADD_TIME
                machine.add_points(z,x,t,tbar,eps,eps_bar,id);
                XPRINT_END_ADD_TIME

                XPRINT_START_OPT_TIME
                FN_RECORD_START
                machine.optimise();
                FN_RECORD_STOP
                XPRINT_END_OPT_TIME

                break;
            }

            case 4:
            {
                long dim;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "Dimension: "; std::cin >> dim; if ( echo ) { std::cout << dim << "\n"; }
                std::cout << "File name (z x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ifstream gs(filename);

                fVECTOR z('r',many);
                f_fVECTOR x('r',many);
                fVECTOR t('r',many);
                fVECTOR tbar('r',many);
                fVECTOR eps('r',many);
                fVECTOR eps_bar('r',many);
                iVECTOR id('r',many);

                for ( i = 1 ; i <= many ; i++ )
                {
                    gs >> z[i-1];

                    (x[i-1]).make_normal(DO_FORCE);
                    (x[i-1]).pad_vector(dim);

                    for ( j = 1 ; j <= dim ; j++ )
                    {
                        gs >> (x[i-1])[j-1];
                    }

                    t[i-1] = 1.0;
                    tbar[i-1] = 1.0;
                    eps[i-1] = 1.0;
                    eps_bar[i-1] = 1.0;

                    id[i-1] = id_count;

                    id_count++;
                }

                gs.close();

                machine.add_points(z,x,t,tbar,eps,eps_bar,id);

                break;
            }

            case 5:
            {
                long dim;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "File name (z x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ifstream gs(filename);

                fVECTOR z('r',many);
                f_fVECTOR x('r',many);
                fVECTOR t('r',many);
                fVECTOR tbar('r',many);
                fVECTOR eps('r',many);
                fVECTOR eps_bar('r',many);
                iVECTOR id('r',many);

                for ( i = 1 ; i <= many ; i++ )
                {
                    gs >> z[i-1];

                    gs >> dim;

                    (x[i-1]).make_normal(DO_FORCE);
                    (x[i-1]).pad_vector(dim);

                    for ( j = 1 ; j <= dim ; j++ )
                    {
                        gs >> (x[i-1])[j-1];
                    }

                    t[i-1] = 1.0;
                    tbar[i-1] = 1.0;
                    eps[i-1] = 1.0;
                    eps_bar[i-1] = 1.0;

                    id[i-1] = id_count;

                    id_count++;
                }

                gs.close();

                XPRINT_START_ADD_TIME
                machine.add_points(z,x,t,tbar,eps,eps_bar,id);
                XPRINT_END_ADD_TIME

                XPRINT_START_OPT_TIME
                FN_RECORD_START
                machine.optimise();
                FN_RECORD_STOP
                XPRINT_END_OPT_TIME

                break;
            }

            case 6:
            {
                long dim;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "File name (z x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ifstream gs(filename);

                fVECTOR z('r',many);
                f_fVECTOR x('r',many);
                fVECTOR t('r',many);
                fVECTOR tbar('r',many);
                fVECTOR eps('r',many);
                fVECTOR eps_bar('r',many);
                iVECTOR id('r',many);

                for ( i = 1 ; i <= many ; i++ )
                {
                    gs >> z[i-1];

                    gs >> dim;

                    (x[i-1]).make_normal(DO_FORCE);
                    (x[i-1]).pad_vector(dim);

                    for ( j = 1 ; j <= dim ; j++ )
                    {
                        gs >> (x[i-1])[j-1];
                    }

                    t[i-1] = 1.0;
                    tbar[i-1] = 1.0;
                    eps[i-1] = 1.0;
                    eps_bar[i-1] = 1.0;

                    id[i-1] = id_count;

                    id_count++;
                }

                gs.close();

                machine.add_points(z,x,t,tbar,eps,eps_bar,id);

                break;
            }

            case 7:
            {
                L_DOUBLE C_new;

                std::cout << "C: "; std::cin >> C_new; if ( echo ) { std::cout << C_new << "\n"; }

                machine.set_C(C_new);
                machine.optimise();

                break;
            }

            case 8:
            {
                L_DOUBLE C_new;

                std::cout << "C+ scale: "; std::cin >> C_new; if ( echo ) { std::cout << C_new << "\n"; }

                machine.set_C_plus_scale(C_new);
                machine.optimise();

                break;
            }

            case 9:
            {
                L_DOUBLE C_new;

                std::cout << "C- scale: "; std::cin >> C_new; if ( echo ) { std::cout << C_new << "\n"; }

                machine.set_C_minus_scale(C_new);
                machine.optimise();

                break;
            }

            case 10:
            {
                L_DOUBLE E_new;

                std::cout << "E: "; std::cin >> E_new; if ( echo ) { std::cout << E_new << "\n"; }

                machine.set_E(E_new);
                machine.optimise();

                break;
            }

            case 11:
            {
                L_DOUBLE E_new;

                std::cout << "E+ scale: "; std::cin >> E_new; if ( echo ) { std::cout << E_new << "\n"; }

                machine.set_E_plus_scale(E_new);
                machine.optimise();

                break;
            }

            case 12:
            {
                L_DOUBLE E_new;

                std::cout << "E- scale: "; std::cin >> E_new; if ( echo ) { std::cout << E_new << "\n"; }

                machine.set_E_minus_scale(E_new);
                machine.optimise();

                break;
            }

            case 13:
            {
                L_DOUBLE nu_new;

                std::cout << "nu: "; std::cin >> nu_new; if ( echo ) { std::cout << nu_new << "\n"; }

                machine.set_nu(nu_new);
                machine.optimise();

                break;
            }

            case 14:
            {
                kernel_dicto K;

                std::cout << "Enter kernel details (sparsity must match training data).\n\n";

                K.prime_io(&(std::cout),echo);

                std::cin >> K;

                if ( (machine.get_kernel()).is_sparse() != K.is_sparse() )
                {
                    std::cout << "Sparsity mismatch - kernel unchanged.\n";
                }

                else
                {
                    machine.set_kernel(K);
                    machine.optimise();
                }

                break;
            }

            case 15:
            {
                std::cout << "Which point: "; std::cin >> i; if ( echo ) { std::cout << i << "\n"; }

                machine.del_point(i);
                machine.optimise();

                break;
            }

            case 16:
            {
                fVECTOR x;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }

                for ( i = 1 ; i <= many ; i++ )
                {
                    if ( echo == 1 ) { x.prime_io(&(std::cout),1); }
                    else             { x.prime_io(&(std::cout),0); }

                    std::cout << "x:\n\n"; std::cin >> x;

                    std::cout << "g(x) = " << machine.test_point(x) << "\n\n";
                }

                break;
            }

            case 17:
            {
                fVECTOR x;
                L_DOUBLE z;
                L_DOUBLE sse;
                L_DOUBLE mse;
                L_DOUBLE rmse;
                L_DOUBLE g;

                sse = 0.0;
                mse = 0.0;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }

                for ( i = 1 ; i <= many ; i++ )
                {
                    std::cout << "z: "; std::cin >> z; if ( echo ) { std::cout << z << "\n"; }

                    if ( echo == 1 ) { x.prime_io(&(std::cout),1); }
                    else             { x.prime_io(&(std::cout),0); }

                    std::cout << "x:\n\n"; std::cin >> x;

                    g = machine.test_point(x);

                    std::cout << "g(x) = " << g << "\n\n";

                    g -= z;
                    g *= g;

                    sse += g;

                    g /= ((L_DOUBLE) many);

                    mse += g;
                }

                rmse = sqrt(mse);

                std::cout << "Sum-squared-error:       " << sse  << "\n";
                std::cout << "Mean-squared-error:      " << mse  << "\n";
                std::cout << "Root-mean-squared-error: " << rmse << "\n\n";

                break;
            }

            case 18:
            {
                long dim;
                long startpoint;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "Number of points to ignore at start: "; std::cin >> startpoint; if ( echo ) { std::cout << startpoint << "\n"; }
                std::cout << "Dimension: "; std::cin >> dim; if ( echo ) { std::cout << dim << "\n"; }
                std::cout << "Input file name (x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ifstream gs(filename);

                std::cout << "Report file name (z format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ofstream hs(filename);

                fVECTOR x('r',dim);

                for ( i = 1 ; i <= many ; i++ )
                {
                    if ( i > startpoint )
                    {
                        for ( j = 1 ; j <= dim ; j++ )
                        {
                            gs >> x[j-1];
                        }
                    }

                    hs << machine.test_point(x) << "\n";
                }

                gs.close();
                hs.close();

                break;
            }

            case 19:
            {
                long dim;
                long startpoint;
                long numpoints;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "Number of points to ignore at start: "; std::cin >> startpoint; if ( echo ) { std::cout << startpoint << "\n"; }
                std::cout << "Dimension: "; std::cin >> dim; if ( echo ) { std::cout << dim << "\n"; }
                std::cout << "Input file name (z x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                numpoints = many-startpoint;

                std::ifstream gs(filename);

                std::cout << "Report file name (z format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ofstream hs(filename);

                fVECTOR x('r',dim);
                L_DOUBLE sse;
                L_DOUBLE mse;
                L_DOUBLE rmse;
                L_DOUBLE z;
                L_DOUBLE g;

                sse = 0;
                mse = 0;

                for ( i = 1 ; i <= many ; i++ )
                {
                    gs >> z;

                    for ( j = 1 ; j <= dim ; j++ )
                    {
                        gs >> x[j-1];
                    }

                    if ( i > startpoint )
                    {
                        g = machine.test_point(x);

                        hs << g << "\n";

                        g -= z;
                        g *= g;

                        sse += g;

                        g /= ((L_DOUBLE) numpoints);

                        mse += g;
                    }
                }

                rmse = sqrt(mse);

                std::cout << "Sum-squared-error:       " << sse  << "\n";
                std::cout << "Mean-squared-error:      " << mse  << "\n";
                std::cout << "Root-mean-squared-error: " << rmse << "\n\n";

                gs.close();
                hs.close();

                break;
            }

            case 20:
            {
                long dim;
                long startpoint;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "Number of points to ignore at start: "; std::cin >> startpoint; if ( echo ) { std::cout << startpoint << "\n"; }
                std::cout << "Input file name (x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ifstream gs(filename);

                std::cout << "Report file name (z format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ofstream hs(filename);

                fVECTOR x('r',dim);

                for ( i = 1 ; i <= many ; i++ )
                {
                    gs >> dim;

                    if ( i > startpoint )
                    {
                        for ( j = 1 ; j <= dim ; j++ )
                        {
                            gs >> x[j-1];
                        }
                    }

                    hs << machine.test_point(x) << "\n";
                }

                gs.close();
                hs.close();

                break;
            }

            case 21:
            {
                long dim;
                long startpoint;
                long numpoints;

                std::cout << "Number of points: "; std::cin >> many; if ( echo ) { std::cout << many << "\n"; }
                std::cout << "Number of points to ignore at start: "; std::cin >> startpoint; if ( echo ) { std::cout << startpoint << "\n"; }
                std::cout << "Input file name (z x format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                numpoints = many-startpoint;

                std::ifstream gs(filename);

                std::cout << "Report file name (g(x) format): "; std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ofstream hs(filename);

                fVECTOR x('r',10);
                L_DOUBLE sse;
                L_DOUBLE mse;
                L_DOUBLE rmse;
                L_DOUBLE z;
                L_DOUBLE g;

                sse = 0;
                mse = 0;

                for ( i = 1 ; i <= many ; i++ )
                {
                    gs >> z;

                    gs >> dim;

                    x.make_zero(DO_FORCE);
                    x.make_normal(DO_FORCE);
                    x.pad_vector(dim);

                    for ( j = 1 ; j <= dim ; j++ )
                    {
                        gs >> x[j-1];
                    }

                    if ( i > startpoint )
                    {
                        g = machine.test_point(x);

                        hs << g << "\n";

                        g -= z;
                        g *= g;

                        sse += g;

                        g /= ((L_DOUBLE) numpoints);

                        mse += g;
                    }
                }

                rmse = sqrt(mse);

                std::cout << "Sum-squared-error:       " << sse  << "\n";
                std::cout << "Mean-squared-error:      " << mse  << "\n";
                std::cout << "Root-mean-squared-error: " << rmse << "\n\n";

                gs.close();
                hs.close();

                break;
            }

            case 22:
            {
                std::cout << machine << "\n";

                break;
            }

            case 23:
            {
                std::cout << "File name: ";
                std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ofstream fs(filename);

                fs << machine;

                fs.close();

                break;
            }

            case 24:
            {
                std::cout << "File name: ";
                std::cin >> filename; if ( echo ) { std::cout << filename << "\n"; }

                std::ifstream gs(filename);

                gs >> machine;

                gs.close();

                break;
            }

            case 25:
            {
                std::cout << "Size:                          " << machine.get_N()   << "\n";
                std::cout << "Support vectors:               " << machine.get_N_S() << "\n";
                std::cout << "Unconstrained support vectors: " << machine.get_N_F() << "\n";

                break;
            }

            case 32:
            {
                XPRINT_START_OPT_TIME
                FN_RECORD_START
                machine.optimise();
                FN_RECORD_STOP
                XPRINT_END_OPT_TIME

                break;
            }

            case 33:
            {
                L_DOUBLE local_tol;

                std::cout << "Tolerance: "; std::cin >> local_tol; if ( echo ) { std::cout << local_tol << "\n"; }

                XPRINT_START_OPT_TIME
                FN_RECORD_START
                machine.optimise(local_tol);
                FN_RECORD_STOP
                XPRINT_END_OPT_TIME

                break;
            }

            case 34:
            {
                L_DOUBLE local_tol;
                int max_iter;

                std::cout << "Tolerance: ";      std::cin >> local_tol; if ( echo ) { std::cout << local_tol << "\n"; }
                std::cout << "Max iterations: "; std::cin >> max_iter;  if ( echo ) { std::cout << max_iter  << "\n"; }

                XPRINT_START_OPT_TIME
                FN_RECORD_START
                machine.optimise_from_zero(local_tol,max_iter);
                FN_RECORD_STOP
                XPRINT_END_OPT_TIME

                break;
            }

            default:
            {
                return;

                break;
            }
        }

        goto repeat;
    }

    return;
}


#ifdef USEALLEGRO
void periodic_report_svc(void *param)
{
    threadtime += *(glob_report_interval);

    std::cerr << threadtime << " sec:  \t";
    std::cerr << "Iterations: " << (*((SVM_pattern *) param)).get_iter();
    std::cerr << "  \tObjective: " << (*((SVM_pattern *) param)).get_J() << "\n";

    return;
}

void periodic_report_svr(void *param)
{
    threadtime += *(glob_report_interval);

    std::cerr << threadtime << " sec:  \t";
    std::cerr << "Iterations: " << (*((SVM_regress *) param)).get_iter();
    std::cerr << "  \tObjective: " << (*((SVM_regress *) param)).get_J() << "\n";

    return;
}
#endif
