
/*
 *  RIMElib: RuntIme Mathematical Equation Library (text processing)
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


#ifndef _mathtext_h
#define _mathtext_h


/*********************************************************************

Mathematical text-processing
============================

Given a null terminated string representing an equation, the function
convert_equation does the following, in order:

1. Removes all whitespace from the equation.
2. Simplify all + and - signs.  Basically, it pairs it down until only
   binary additions and unary negations remain.
3. Convert all numbers n into the form one[n].
4. Standardise bracketting.
5. Eliminate unary negations.
6. Convert a^b (a to the power of b) to pow[1](a,b), right to left.
7. Convert a*b and a/b to mul[1](a,b) and div[1](a,b),
   respectively, working from left to right.
8. Convert a+b to add[1](a,b) working from left to right.

and sets result to point to a new string (the old remains unchanged)
in a form ready for processing.  It will return either the length of
the result (>= 0) or a negative error code if problems are encountered
(see defined syntax errors of form MERR_SYNTAX_ERROR_... below).  The
internal functions used are, in order:

1. remove_whitespace
2. simplify_signs
3. convert_numbers
4. convert_bracks
5. convert_standbrack
6. convert_negations
7. convert_pow
8. convert_muldiv
9. convert_add
10. convert_passes

NB: - neither e nor E may be used as a function on its own.
    - the function pass is reserved for internal use.


Internal functions
==================

pre-processing: The following functions modify the string express of length
                len, returning the new length.

remove_whitespace - removes whitespace from express
simplify_signs    - remove sign conglomerates, remove redundant unary
                    posations and convert binary subtractions to unary
                    negations followed by addition.  This function assumes
                    all whitespace has been removed before calling.


analysis: The following functions analyse the string express of length len
          to obtain information.  They assume that no whitespace is present
          in the string, that all +'s are binary additions and that all -'s
          are unary negations.

len_to_nonalpha - find number of chars before next non alpha character
                  return 0 if there are no alpha's before the 1st non-alhpa
                  return len if no non-alpha's are found
                  note that _ is considered an alpha character
len_to_open     - find number of chars before next (,[ or {
                  return 0 if string starts with said (obviously)
                  return len if none of said present
len_to_br_open  - find number of chars before next (
                  return 0 if string starts with said (obviously)
                  return len if none of said present
len_to_br_close - find number of chars before next ) at current level,
                  do basic syntax checking
                  error if no ) found before end
len_to_sq_close - find number of chars before next ] at current level,
                  do basic syntax checking
                  error if no ] found before end
len_to_comma    - find number of chars before next , at current level
len_to_comma_br - find number of chars before next , or ) at current level
len_to_uint_end - find number of chars in unsigned int (0 converted to error)
len_to_int_end  - find number of chars in signed int (0 converted to error)
len_to_real_end - find number of chars in real (0 converted to error)


further processing: same assumptions, create new strings.  Each converter
                    assumes previous converters have already been called.

convert_numbers    - takes a string and creates a new (null-terminated)
                     string where all numbers n are replaced with the
                     expression one[n] (this can be changed using the macros
                     MATHS_CONST_LEFT and MATHS_CONST_RIGHT).  Actual
                     operation is to process from left to right thusly:
                     - on finding a [ skip to ].
                     - on finding a *+/() or alpha go on to the next char.
                     - on finding a letter skip to the end of the function.
                     - on finding a 0123456789-. then read the number (a
                       number may include 0123456789-.eE but may not start
                       with eE) and enclose it appropriately.
convert_bracks     - convert "bare" brackets onto proper expessions.  A bare
                     bracket can take any of the following forms:
                     - ( at the start of the function
                     - ( preceeded by a non-alpha character not including ]})
convert_standbrack - standardise brackets.  All functions should be a series
                     of alphas followed by [...](...) in that order.
                     This function assumes the order is correct already, but
                     that some of the bracket pairs may be missing.
convert_negations  - remove unary negations by converting putting them inside
                     multiplier brackers [].
convert_pow        - replace a^b with pow[1](a,b) from right to left
convert_muldiv     - replace a*b with mul[1](a,b) from left to right and
                     a/b with mul[1](a,inv(b)) left to right simultaneously
convert_add        - replace a+b with add[1](a,b) from left to right
convert_passes     - replace expressions pass[1](...) with ...

Generic subroutines:

convert_ltor(int lev_size, char *sep_expr, char **start_expr, char **mid_expr, char **end_expr, char *what, int *len);
convert_rtol(int lev_size, char *sep_expr, char **start_expr, char **mid_expr, char **end_expr, char *what, int *len);

*********************************************************************/

long convert_equation(char *express, char **result);

long remove_whitespace(char **express, long len);
long simplify_signs(char **express, long len);

long len_to_nonalpha(char *express, long len);
long len_to_open(char *express, long len);
long len_to_br_open(char *express, long len);
long len_to_br_close(char *express, long len);
long len_to_sq_close(char *express, long len);
long len_to_comma(char *express, long len);
long len_to_comma_br(char *express, long len);
long len_to_uint_end(char *express, long len);
long len_to_int_end(char *express, long len);
long len_to_real_end(char *express, long len);

long convert_numbers(char **express, long len);
long convert_bracks(char **express, long len);
long convert_standbrack(char **express, long len);
long convert_negations(char **express, long len);

long convert_pow(char **express, long len);
long convert_muldiv(char **express, long len);
long convert_add(char **express, long len);
long convert_passes(char **express, long len);

char *convert_ltor(int lev_size, char *sep_expr, char **start_expr, char **mid_expr, char **end_expr, char *what, long *len);
char *convert_rtol(int lev_size, char *sep_expr, char **start_expr, char **mid_expr, char **end_expr, char *what, long *len);


/*

Macros

*/

#define MATHS_CONST_LEFT        "one["
#define MATHS_CONST_RIGHT       "]()"

#define MATHS_CONST_BARE_BRACK  "pass"

#define MATHS_POW_STRING_LEFT   "pow[1]("
#define MATHS_POW_STRING_MID    ","
#define MATHS_POW_STRING_RIGHT  ")"

#define MATHS_MUL_STRING_LEFT   "mul[1]("
#define MATHS_MUL_STRING_MID    ","
#define MATHS_MUL_STRING_RIGHT  ")"

#define MATHS_DIV_STRING_LEFT   "div[1]("
#define MATHS_DIV_STRING_MID    ","
#define MATHS_DIV_STRING_RIGHT  ")"

#define MATHS_ADD_STRING_LEFT   "add[1]("
#define MATHS_ADD_STRING_MID    ","
#define MATHS_ADD_STRING_RIGHT  ")"



/*

Syntax errors:

*/

#define MERR_SYNTAX_ERROR_SIGN_AT_END_A         -1
#define MERR_SYNTAX_ERROR_SIGN_AT_END_B         -2
#define MERR_SYNTAX_ERROR_E_AT_START_A          -3
#define MERR_SYNTAX_ERROR_E_AT_START_B          -4
#define MERR_SYNTAX_ERROR_NO_ALPHAS             -5
#define MERR_SYNTAX_ERROR_UNPAIRED_BRACKETS_A   -6
#define MERR_SYNTAX_ERROR_UNPAIRED_BRACKETS_B   -7
#define MERR_SYNTAX_ERROR_UNPAIRED_BRACKETS_C   -8
#define MERR_SYNTAX_ERROR_UNPAIRED_BRACKETS_D   -9
#define MERR_SYNTAX_ERROR_UNPAIRED_BRACKETS_E   -10
#define MERR_SYNTAX_ERROR_UNPAIRED_BRACKETS_F   -11
#define MERR_SYNTAX_ERROR_COMMA_NOT_FOUND       -12
#define MERR_SYNTAX_ERROR_ZERO_LENGTH_UINT      -13
#define MERR_SYNTAX_ERROR_ZERO_LENGTH_INT       -14
#define MERR_SYNTAX_ERROR_ZERO_LENGTH_REAL      -15
#define MERR_SYNTAX_ERROR_MISPLACED_NEGATION    -16
#define MERR_SYNTAX_ERROR_E_SINGLETON_A         -17
#define MERR_SYNTAX_ERROR_E_SINGLETON_B         -18
#define MERR_SYNTAX_ERROR_E_SINGLETON_C         -19
#define MERR_SYNTAX_ERROR_E_SINGLETON_D         -20
#define MERR_SYNTAX_ERROR_NEG_SINGLETON_A       -21
#define MERR_SYNTAX_ERROR_NEG_SINGLETON_B       -22
#define MERR_SYNTAX_ERROR_NEG_SINGLETON_C       -23
#define MERR_SYNTAX_ERROR_NEG_SINGLETON_D       -24
#define MERR_SYNTAX_ERROR_REAL_DOT_SINGLETON    -25
#define MERR_SYNTAX_ERROR_MISSED_NUMBER_A       -26
#define MERR_SYNTAX_ERROR_MISSED_NUMBER_B       -27
#define MERR_SYNTAX_ERROR_DOUBLE_DOT            -28
#define MERR_SYNTAX_ERROR_REAL_MISPLACED_OBJECT -29
#define MERR_SYNTAX_ERROR_E_END_SINGLETON       -30
#define MERR_SYNTAX_ERROR_NOTHING               -31
#define MERR_SYNTAX_ERROR_NO_FN_A               -32
#define MERR_SYNTAX_ERROR_NO_FN_B               -33
#define MERR_SYNTAX_ERROR_SHORT_FN              -34
#define MERR_SYNTAX_ERROR_NEG_END               -35
#define MERR_MEM_ERROR_A                        -36
#define MERR_MEM_ERROR_B                        -37
#define MERR_MEM_ERROR_C                        -38
#define MERR_MEM_ERROR_D                        -39
#define MERR_MEM_ERROR_E                        -40
#define MERR_MEM_ERROR_F                        -41
#define MERR_MEM_ERROR_G                        -42
#define MERR_MEM_ERROR_H                        -43
#define MERR_MEM_ERROR_I                        -44
#define MERR_MEM_ERROR_J                        -45
#define MERR_MEM_ERROR_K                        -46
#define MERR_MEM_ERROR_L                        -47
#define MERR_MEM_ERROR_M                        -48
#define MERR_MEM_ERROR_N                        -49
#define MERR_MEM_ERROR_O                        -50
#define MERR_MEM_ERROR_P                        -51
#define MERR_MEM_ERROR_Q                        -52
#define MERR_MEM_ERROR_R                        -53
#define MERR_MEM_ERROR_S                        -54
#define MERR_MEM_ERRORB_A                       -55
#define MERR_MEM_ERRORB_B                       -56
#define MERR_MEM_ERRORB_C                       -57
#define MERR_MEM_ERRORB_D                       -58
#define MERR_MEM_ERRORB_E                       -59
#define MERR_MEM_ERRORB_F                       -60
#define MERR_MEM_ERRORB_G                       -61
#define MERR_MEM_ERRORB_H                       -62
#define MERR_MEM_ERRORB_I                       -63
#define MERR_MEM_ERRORB_J                       -64
#define MERR_MEM_ERRORB_K                       -65
#define MERR_MEM_ERRORB_L                       -66
#define MERR_MEM_ERRORB_M                       -67
#define MERR_SYNTAX_ERROR_PASSREM_A             -68
#define MERR_SYNTAX_ERROR_PASSREM_B             -69
#define MERR_SYNTAX_ERROR_PASSREM_C             -70
#define MERR_SYNTAX_ERROR_PASSREM_D             -71
#define MERR_SYNTAX_ERROR_PASSREM_E             -72
#define MERR_SYNTAX_ERROR_PASSREM_F             -73
#define MERR_SYNTAX_ERROR_PASSREM_G             -74
#define MERR_SYNTAX_ERROR_PASSREM_H             -75
#define MERR_SYNTAX_ERROR_PASSREM_I             -76
#define MERR_SYNTAX_ERROR_PASSREM_J             -77
#define MERR_SYNTAX_ERROR_PASSREM_K             -78
#define MERR_SYNTAX_ERROR_PASSREM_L             -79
#define MERR_SYNTAX_ERROR_PASSREM_M             -80
#define MERR_SYNTAX_ERROR_PASSREM_N             -81
#define MERR_SYNTAX_ERROR_PASSREM_O             -82
#define MERR_SYNTAX_ERROR_BADCHAR               -83
#define MERR_CONCLUDING_MINUS                   -84


#define MERR_MAX_ERROR                          -90



#endif

