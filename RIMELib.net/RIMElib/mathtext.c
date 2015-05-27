#include "stdafx.h"
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


#include <string.h>
#include <malloc.h>
#include <stdlib.h>

#include "mathtext.h"


void mathtext_exitclean(void);



/*********************************************************************

Mathematical text-processing
============================

*********************************************************************/

/*

convert_equation - this is a shortcut - it re-allocates the string and calls
                   the following, in order: remove_whitespace,
                   simplify_signs, convert_numbers, convert_bracks,
                   convert_standbrack, convert_negations, convert_pow,
                   convert_muldiv and convert_add.

*/

long convert_equation(char *express, char **result)
{
    char *xresult[1];
    long len;

    if ( ( xresult[0] = (char *) malloc((strlen(express)+1)*sizeof(char)) ) == NULL )
    {
        return MERR_MEM_ERROR_A;
    }

    strcpy(xresult[0],express);

    len = strlen(xresult[0]);

    if ( len > 0 ) { len = remove_whitespace(xresult,len);      }
    if ( len > 0 ) { len = simplify_signs(xresult,len);         }
    if ( len > 0 ) { len = convert_numbers(xresult,len);        }
    if ( len > 0 ) { len = convert_bracks(xresult,len);         }
    if ( len > 0 ) { len = convert_standbrack(xresult,len);     }
    if ( len > 0 ) { len = convert_negations(xresult,len);      }
    if ( len > 0 ) { len = convert_pow(xresult,len);            }
    if ( len > 0 ) { len = convert_muldiv(xresult,len);         }
    if ( len > 0 ) { len = convert_add(xresult,len);            }
    if ( len > 0 ) { len = convert_passes(xresult,len);         }

    if ( len > 0 ) { xresult[0][len] = '\0'; }

    *result = xresult[0];

    return len;
}




/*

remove_whitespace - removes whitespace from express.  Also checks to
                    make sure that all characters are legal.

*/

long remove_whitespace(char **express, long len)
{
    long i = 0;
    long j;

    while ( i < len )
    {
        switch ( (*express)[i] )
        {
            case ' ':
            case '\t':
            {
                j = i+1;

                while ( j < len )
                {
                    (*express)[j-1] = (*express)[j];

                    j++;
                }

                i--;

                len--;

                break;
            }

            case '0': case '1': case '2': case '3': case '4':
            case '5': case '6': case '7': case '8': case '9':
            case ',': case '.': case '/': case '[': case ']':
            case '^': case '*': case '(': case ')': case '-':
            case '+': case '_':
            case 'a': case 'b': case 'c': case 'd': case 'e':
            case 'f': case 'g': case 'h': case 'i': case 'j':
            case 'k': case 'l': case 'm': case 'n': case 'o':
            case 'p': case 'q': case 'r': case 's': case 't': 
            case 'u': case 'v': case 'w': case 'x': case 'y': 
            case 'z':
            case 'A': case 'B': case 'C': case 'D': case 'E':
            case 'F': case 'G': case 'H': case 'I': case 'J':
            case 'K': case 'L': case 'M': case 'N': case 'O':
            case 'P': case 'Q': case 'R': case 'S': case 'T': 
            case 'U': case 'V': case 'W': case 'X': case 'Y': 
            case 'Z':
            {
                break;
            }

            default:
            {
                return MERR_SYNTAX_ERROR_BADCHAR;

                break;
            }
        }

        i++;
    }

    (*express)[len] = '\0';

    return len;
}




/*

simplify_signs    - remove sign conglomerates, remove redundant unary
                    posations and convert binary subtractions to unary
                    negations followed by addition.  This function assumes
                    all whitespace has been removed before calling.

Method:

This is a three pass function.

Pass 1: On the first pass, sign conglomerates (for example +---+++-, or more
        likely +-, are removed, replaced by either + (if the number of -'s
        in the conglomerate is even) or - (if the number of -'s is odd).
        Also checks to ensure that no expression (or sub expression) ends in
        a stray + or -.
        This pass operates searching for pairs ++, -- (both replaced by +)
        and +-, -+ (both replaced by -).  When a pair is found, it is cut
        out and replaced, and the search recommences on the same character
        (to find string longer than 2).

Pass 2: Removes all redundant posations.  For example, sin(+x) is replaced
        by sin(x).  This simplifies the job (done later) of replacing all
        operators by functions.  Rather than having two operations using +
        (binary addition and unary posation), removing all posations (which
        are redundant in any case) means that any + can be assumed to be a
        binary operator.
        A redundant posation can be seen to be so if it follows a bracket
        ( or [, an unambiguous binary operator , * / ^ or and e or E which
        is preceeded by a digit 0,1,2,3,4,5,6,7,8,9 or a . to indicate that
        this e or E is a exponentiation indicator.  When a posation is found
        it is simply removed.

Pass 3: Replace binary subtraction with unary negation of the second element
        followed by binary addition - that is, replace - with +-.
        When a - is encountered, it is first checked to see if it is a unary
        negation (using the same test as for unary posation).  If this test
        fails (ie. it is a binary subtraction operator) then a + is inserted
        into the string before the -.
        This operation simplifies the later job of replacing all operators
        with functions.  Once this is completed (along with pass 2) then it
        is possible to simply *assume* that all +'s are binary addition
        operators, and all -'s are unary subtraction operators, thus
        removing unnecessary complications later.

*/

long simplify_signs(char **express, long len)
{
    long i;
    long j;
    char *result;

    /* First pass - remove conglomerates (eg. +--+++-). */

    i = 0;

    while ( i < len )
    {
        switch ( (*express)[i] )
        {
            case '+':
            {
                if ( i >= len-1 )
                {
                    /*
                        Cannot have a + or - at the end of an equation -
                        this is bad syntax.  For eg.  sin(x)+ is clearly
                        meaningless
                    */

                    return MERR_SYNTAX_ERROR_SIGN_AT_END_A;
                }

                switch ( (*express)[i+1] )
                {
                    case '+':
                    case '-':
                    {
                        /*
                           +- is replaced by -
                           ++ is replaced by +

                           i is decremented here so that this char will be
                           checked again.  Thus arbitrarily long sequences
                           of + and - can be dealt with.
                        */

                        j = i+1;

                        while ( j < len )
                        {
                            (*express)[j-1] = (*express)[j];

                            j++;
                        }

                        i--;

                        len--;

                        break;
                    }

                    default:
                    {
                        break;
                    }
                }

                break;
            }

            case '-':
            {
                if ( i >= len-1 )
                {
                    /*
                        Cannot have a + or - at the end of an equation -
                        this is bad syntax.  For eg.  sin(x)+ is clearly
                        meaningless
                    */

                    return MERR_SYNTAX_ERROR_SIGN_AT_END_B;
                }

                switch ( (*express)[i+1] )
                {
                    case '+':
                    {
                        /*
                           -+ is replaced by -

                           i is decremented here so that this char will be
                           checked again.  Thus arbitrarily long sequences
                           of + and - can be dealt with.
                        */

                        j = i+1;

                        (*express)[j] = '-';

                        while ( j < len )
                        {
                            (*express)[j-1] = (*express)[j];

                            j++;
                        }

                        i--;

                        len--;

                        break;
                    }

                    case '-':
                    {
                        /*
                           -- is replaced by +

                           i is decremented here so that this char will be
                           checked again.  Thus arbitrarily long sequences
                           of + and - can be dealt with.
                        */

                        j = i+1;

                        (*express)[j] = '+';

                        while ( j < len )
                        {
                            (*express)[j-1] = (*express)[j];

                            j++;
                        }

                        i--;

                        len--;

                        break;
                    }

                    default:
                    {
                        break;
                    }
                }

                break;
            }

            default:
            {
                break;
            }
        }

        i++;
    }

    /* Second pass - remove redundant posations (unary + is redundant) */

    i = 0;

    while ( i < len )
    {
        switch ( (*express)[i] )
        {
            case '+':
            {
                j = 0;

                /*
                   Set j if this is a redundant posation.
                */

                if ( i )
                {
                    switch ( (*express)[i-1] )
                    {
                        case '(':
                        case '[':
                        case ',':
                        case '*':
                        case '/':
                        case '^':
                        {
                            /*
                               If + follows any of these it can be assumed
                               to be a redundant posation.
                            */

                            j = 1;

                            break;
                        }

                        case 'e':
                        case 'E':
                        {
                            /*
                               If + was preceeded by e or E then it may be
                               a redundant posation if e was acting as an
                               exponent operation in a real number.
                            */

                            if ( i <= 1 )
                            {
                                /*
                                   e and E are on there own cannot be
                                   used as a function.
                                */

                                return MERR_SYNTAX_ERROR_E_AT_START_A;
                            }

                            /*
                               Test to see what came before the e or E to
                               see if this is a number or not.  If it is a
                               number then + must be a redundant posation.
                            */

                            switch ( (*express)[i-2] )
                            {
                                case '0':
                                case '1':
                                case '2':
                                case '3':
                                case '4':
                                case '5':
                                case '6':
                                case '7':
                                case '8':
                                case '9':
                                case '.':
                                {
                                    j = 1;

                                    break;
                                }

                                default:
                                {
                                    break;
                                }
                            }

                            break;
                        }

                        default:
                        {
                            break;
                        }
                    }
                }

                else
                {
                    /*
                       A + at the start of an expression can be assumed to
                       be a redundant posation.
                    */

                    j = 1;
                }

                if ( j )
                {
                    /*
                       Remove the redundant posation.  As conglomerates have
                       already been removed, we can assume that the next
                       character is not a +.
                    */

                    j = i+1;

                    while ( j < len )
                    {
                        (*express)[j-1] = (*express)[j];

                        j++;
                    }

                    len--;
                }

                break;
            }

            default:
            {
                break;
            }
        }

        i++;
    }

    /* Third pass - replace binary subtraction with +- sequence (we do  */
    /*              this to make parsing easier.  Having done this, we  */
    /*              may assume that all +'s are binary additions, and   */
    /*              that all -'s are unary negations).                  */

    i = 0;

    while ( i < len )
    {
        switch ( (*express)[i] )
        {
            case '-':
            {
                j = 1;

                /*
                   Reset j if this is not a binary negation.
                */

                if ( i )
                {
                    switch ( (*express)[i-1] )
                    {
                        case '(':
                        case '[':
                        case ',':
                        case '*':
                        case '/':
                        case '^':
                        {
                            /*
                               If the - follows any of these then it must
                               be a unary negation, so we can ignore it.
                            */

                            j = 0;

                            break;
                        }

                        case 'e':
                        case 'E':
                        {
                            if ( i <= 1 )
                            {
                                /*
                                   e and E are on there own cannot be
                                   used as a function.
                                */

                                return MERR_SYNTAX_ERROR_E_AT_START_B;
                            }

                            /*
                               As when removing redundant posations, if -
                               follows an e or E then we may be dealing with
                               a unary negation *if* the e or E is part of a
                               number.  So test the character before the e
                               or E.
                            */

                            switch ( (*express)[i-2] )
                            {
                                case '0':
                                case '1':
                                case '2':
                                case '3':
                                case '4':
                                case '5':
                                case '6':
                                case '7':
                                case '8':
                                case '9':
                                case '.':
                                {
                                    j = 0;

                                    break;
                                }

                                default:
                                {
                                    break;
                                }
                            }

                            break;
                        }

                        default:
                        {
                            break;
                        }
                    }
                }

                else
                {
                    j = 0;
                }

                if ( j )
                {
                    /*
                       Having found a binary subtraction, replace it with
                       +- (unary negation followed by binary addition).
                    */

                    if ( ( result = (char *) malloc((len+2)*sizeof(char)) ) == NULL )
                    {
                        return MERR_MEM_ERROR_B;
                    }

                    j = 0;

                    while ( j < i )
                    {
                        result[j] = (*express)[j];

                        j++;
                    }

                    result[j] = '+';

                    while ( j < len )
                    {
                        result[j+1] = (*express)[j];

                        j++;
                    }

                    free(*express);

                    *express = result;

                    len++;

                    i++;
                }

                break;
            }

            default:
            {
                break;
            }
        }

        i++;
    }

    (*express)[len] = '\0';

    return len;
}




/*

len_to_nonalpha - find number of chars before next non alpha character.
                  return 0 if there are no alpha's before the 1st non-alhpa
                  return len if no non-alpha's are found.
                  note that _ is considered an alpha character.

*/

long len_to_nonalpha(char *express, long len)
{
    long i = 0;

    while ( i < len )
    {
        switch ( express[i] )
        {
            case 'a': case 'b': case 'c': case 'd': case 'e': case 'f':
            case 'g': case 'h': case 'i': case 'j': case 'k': case 'l':
            case 'm': case 'n': case 'o': case 'p': case 'q': case 'r':
            case 's': case 't': case 'u': case 'v': case 'w': case 'x':
            case 'y': case 'z':
            case 'A': case 'B': case 'C': case 'D': case 'E': case 'F':
            case 'G': case 'H': case 'I': case 'J': case 'K': case 'L':
            case 'M': case 'N': case 'O': case 'P': case 'Q': case 'R':
            case 'S': case 'T': case 'U': case 'V': case 'W': case 'X':
            case 'Y': case 'Z':
            case '_':
            {
                break;
            }

            default:
            {
                return i;

                break;
            }
        }

        i++;
    }

    return len;
}




/*

len_to_open     - find number of chars before next (,[ or {
                  return 0 if string starts with said (obviously)
                  return len if none of said present

*/

long len_to_open(char *express, long len)
{
    long i = 0;

    while ( i < len )
    {
        switch ( express[i] )
        {
            case '(':
            case '[':
            {
                return i;

                break;
            }

            default:
            {
                break;
            }
        }

        i++;
    }

    return len;
}




/*

len_to_br_open  - find number of chars before next (
                  return 0 if string starts with said (obviously)
                  return len if none of said present

*/

long len_to_br_open(char *express, long len)
{
    long i = 0;

    while ( i < len )
    {
        switch ( express[i] )
        {
            case '(':
            {
                return i;

                break;
            }

            default:
            {
                break;
            }
        }

        i++;
    }

    return len;
}




/*

len_to_br_close - find number of chars before next ) at current level,
                  do basic syntax checking
                  error if no ) found before end

*/

long len_to_br_close(char *express, long len)
{
    long i = 0;
    long j;

    while ( i < len )
    {
        switch ( express[i] )
        {
            case '(':
            {
                /*
                   We can pass over whatever is in these brackets - it is
                   not at the current level.
                */

                i++;

                j = len_to_br_close(express+i,len-i);

                if ( j < 0 ) { return j; }

                i += j;

                break;
            }

            case '[':
            {
                /*
                   We can pass over whatever is in these brackets - it is
                   not at the current level.
                */

                i++;

                j = len_to_sq_close(express+i,len-i);

                if ( j < 0 ) { return j; }

                i += j;

                break;
            }

            case ')':
            {
                /*
                   OK, this is it.
                */

                return i;

                break;
            }

            case ']':
            {
                /*
                   As this doesn't line up with any [, we can assume that
                   something is wrong.
                */

                return MERR_SYNTAX_ERROR_UNPAIRED_BRACKETS_A;

                break;
            }

            default:
            {
                break;
            }
        }

        i++;
    }

    return MERR_SYNTAX_ERROR_UNPAIRED_BRACKETS_B;
}




/*

len_to_sq_close - find number of chars before next ] at current level,
                  do basic syntax checking
                  error if no ] found before end

*/

long len_to_sq_close(char *express, long len)
{
    long i = 0;
    long j;

    while ( i < len )
    {
        switch ( express[i] )
        {
            case '(':
            {
                i++;

                j = len_to_br_close(express+i,len-i);

                if ( j < 0 ) { return j; }

                i += j;

                break;
            }

            case '[':
            {
                i++;

                j = len_to_sq_close(express+i,len-i);

                if ( j < 0 ) { return j; }

                i += j;

                break;
            }

            case ')':
            {
                return MERR_SYNTAX_ERROR_UNPAIRED_BRACKETS_C;

                break;
            }

            case ']':
            {
                return i;

                break;
            }

            default:
            {
                break;
            }
        }

        i++;
    }

    return MERR_SYNTAX_ERROR_UNPAIRED_BRACKETS_D;
}




/*

len_to_comma    - find number of chars before next , at current level

*/

long len_to_comma(char *express, long len)
{
    long i = 0;
    long j;

    while ( i < len )
    {
        switch ( express[i] )
        {
            case '(':
            {
                i++;

                j = len_to_br_close(express+i,len-i);

                if ( j < 0 ) { return j; }

                i += j;

                break;
            }

            case '[':
            {
                i++;

                j = len_to_sq_close(express+i,len-i);

                if ( j < 0 ) { return j; }

                i += j;

                break;
            }

            case ')':
            {
                return MERR_SYNTAX_ERROR_UNPAIRED_BRACKETS_E;

                break;
            }

            case ']':
            {
                return MERR_SYNTAX_ERROR_UNPAIRED_BRACKETS_F;

                break;
            }

            case ',':
            {
                return i;

                break;
            }

            default:
            {
                break;
            }
        }

        i++;
    }

    return MERR_SYNTAX_ERROR_COMMA_NOT_FOUND;
}




/*

len_to_comma_br - find number of chars before next , or ) at current level

*/

long len_to_comma_br(char *express, long len)
{
    long i = 0;
    long j;

    while ( i < len )
    {
        switch ( express[i] )
        {
            case '(':
            {
                i++;

                j = len_to_br_close(express+i,len-i);

                if ( j < 0 ) { return j; }

                i += j;

                break;
            }

            case '[':
            {
                i++;

                j = len_to_sq_close(express+i,len-i);

                if ( j < 0 ) { return j; }

                i += j;

                break;
            }

            case ']':
            {
                return MERR_SYNTAX_ERROR_UNPAIRED_BRACKETS_F;

                break;
            }

            case ')':
            case ',':
            {
                return i;

                break;
            }

            default:
            {
                break;
            }
        }

        i++;
    }

    return MERR_SYNTAX_ERROR_COMMA_NOT_FOUND;
}




/*

len_to_uint_end - find number of chars in unsigned int (0 converted to error)

*/

long len_to_uint_end(char *express, long len)
{
    long i = 0;

    while ( i < len )
    {
        switch ( express[i] )
        {
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
            {
                break;
            }

            default:
            {
                if ( i ) { return i; }

                return MERR_SYNTAX_ERROR_ZERO_LENGTH_UINT;

                break;
            }
        }

        i++;
    }

    if ( len ) { return len; }

    return MERR_SYNTAX_ERROR_ZERO_LENGTH_UINT;
}




/*

len_to_int_end  - find number of chars in signed int (0 converted to error)

*/

long len_to_int_end(char *express, long len)
{
    long i;

    if ( len )
    {
        switch ( express[0] )
        {
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
            {
                return len_to_uint_end(express,len);

                break;
            }

            case '-':
            {
                i = len_to_uint_end(express+1,len-1);

                if ( i ) { return i+1; }

                /*
                   A single -, on its own, is not an integer.
                */

                if ( i == MERR_SYNTAX_ERROR_ZERO_LENGTH_UINT )
                {
                    return MERR_SYNTAX_ERROR_ZERO_LENGTH_INT;
                }

                /*
                   Other errors can pass through.
                */

                return i;

                break;
            }

            default:
            {
                return MERR_SYNTAX_ERROR_ZERO_LENGTH_INT;

                break;
            }
        }
    }

    return MERR_SYNTAX_ERROR_ZERO_LENGTH_INT;
}




/*

len_to_real_end - find number of chars in real (0 converted to error)

*/

long len_to_real_end(char *express, long len)
{
    long i = 0;
    long j;
    int no_number = 1;

    /*
       Note: this is reliant on the fact that all -'s are unary negations,
       and hence part of the real number, and that all +'s are binary
       additions, and hence NOT part of the real number.
    */

    if ( len == 0 )
    {
        return MERR_SYNTAX_ERROR_ZERO_LENGTH_REAL;
    }

    /*
       Start by looking at the first letter.
    */

    switch ( express[i] )
    {
        case '-':
        {
            /*
               A real number can start with a - so long as this is not the
               entirety of the number.  Skip over this.
            */

            i++;

            if ( i >= len )
            {
                return MERR_SYNTAX_ERROR_MISPLACED_NEGATION;
            }

            break;
        }

        case 'e':
        case 'E':
        {
            /*
               A real number cannot start with an e or E.
            */

            return MERR_SYNTAX_ERROR_E_SINGLETON_A;

            break;
        }

        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
        case '.':
        default:
        {
            break;
        }
    }

    /*
       Now we may have an unsigned integer - we can move past this.
    */

    switch ( express[i] )
    {
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
        {
            j = len_to_uint_end(express+i,len-i);

            if ( j < 0 ) { return j; }

            i += j;

            if ( i >= len ) { return i; }

            /*
               We need to record if there is a number here.  This will
               effect the syntax of what follows.
            */

            no_number = 0;

            break;
        }

        case 'e':
        case 'E':
        {
            /*
               A number cannot start with e,E,-e or -E.
            */

            return MERR_SYNTAX_ERROR_E_SINGLETON_B;

            break;
        }

        case '-':
        {
            /*
               We already dealt with the possibility of a - at the start
               of the number.  If we found another then there are syntax
               errors here.
            */

            return MERR_SYNTAX_ERROR_NEG_SINGLETON_A;

            break;
        }

        case '.':
        {
            break;
        }

        default:
        {
            if ( express[0] != '-' )
            {
                return MERR_SYNTAX_ERROR_NEG_SINGLETON_B;
            }

            return MERR_SYNTAX_ERROR_ZERO_LENGTH_REAL;

            break;
        }
    }

    /*
       Then we might have a decimal point, followed by an unsigned int.
    */

    switch ( express[i] )
    {
        case '.':
        {
            i++;

            if ( i >= len )
            {
                /*
                   A number ending in a dot is only valid if a number
                   preceeded this.
                */

                if ( no_number )
                {
                    return MERR_SYNTAX_ERROR_REAL_DOT_SINGLETON;
                }

                return i;
            }

            switch ( express[i] )
            {
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                {
                    /*
                       A number following the can be parsed over ok.
                    */

                    j = len_to_uint_end(express+i,len-i);

                    if ( j < 0 ) { return j; }

                    i += j;

                    if ( i >= len ) { return i; }

                    /*
                       Record that a number has happenned here.
                    */

                    no_number = 0;

                    break;
                }

                case '-':
                {
                    /*
                       .- inside a string is not possible.
                    */

                    return MERR_SYNTAX_ERROR_NEG_SINGLETON_C;

                    break;
                }

                case '.':
                {
                    /*
                       .. inside a string is not possible.
                    */

                    return MERR_SYNTAX_ERROR_DOUBLE_DOT;

                    break;
                }

                case 'e':
                case 'E':
                {
                    /*
                       e or E is OK so long as it follows a number (in
                       which case it represents an exponent).  Otherwise,
                       it's syntactically incorrect.
                    */

                    if ( no_number )
                    {
                        return MERR_SYNTAX_ERROR_E_SINGLETON_C;
                    }

                    break;
                }

                default:
                {
                    break;
                }
            }

            break;
        }

        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
        {
            /*
               If we are finding a number here, something has gone very
               wrong, as we just skipped past all numbers (ie. this is
               not possible).
            */

            return MERR_SYNTAX_ERROR_MISSED_NUMBER_A;

            break;
        }

        case '-':
        {
            /*
               Again, if we find this, we have problems.
            */

            return MERR_SYNTAX_ERROR_NEG_SINGLETON_D;

            break;
        }

        case 'e':
        case 'E':
        {
            /*
               e or E is OK so long as it follows a number (in which case
               it represents an exponent).  Otherwise, it's syntactically
               incorrect.
            */

            if ( no_number )
            {
                return MERR_SYNTAX_ERROR_E_SINGLETON_D;
            }

            break;
        }

        default:
        {
            break;
        }
    }

    if ( no_number )
    {
        /*
           If we still haven't found a number then this *isn't* a number,
           so something is wrong.
        */

        return MERR_SYNTAX_ERROR_MISSED_NUMBER_B;
    }

    /* Finally, allow for possible exponent */

    switch ( express[i] )
    {
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
        case '.':
        case '-':
        {
            /*
               Cannot find any of these unless something is seriously wrong.
            */

            return MERR_SYNTAX_ERROR_REAL_MISPLACED_OBJECT;

            break;
        }

        case 'e':
        case 'E':
        {
            i++;

            if ( i >= len )
            {
                /*
                   An exponent indicator must be followed by a number, or
                   this is syntactically invalid.
                */

                return MERR_SYNTAX_ERROR_E_END_SINGLETON;
            }

            j = len_to_real_end(express+i,len-i);

            if ( j < 0 ) { return j; }

            i += j;

            break;
        }

        default:
        {
            break;
        }
    }

    /*
       And we're done!
    */

    return i;
}




/*

convert_numbers - takes a string and creates a new (null-terminated) string
                  where all numbers n are replaced with the expression
                  one[n]() (this can be changed using the macros
                  MATHS_CONST_LEFT and MATHS_CONST_RIGHT).  Actual operation
                  is to process from left to right thusly:
                  - on finding a [ skip to ].
                  - on finding a *+/() or alpha go on to the next character.
                  - on finding a letter skip to the end of the function.
                  - on finding a 0123456789-. then read the number (a number
                    may include 0123456789-.eE but may not start with eE)
                    and enclose it appropriately.

*/

long convert_numbers(char **express, long len)
{
    long i = 0;
    long j;
    long k;
    char left_string[] = MATHS_CONST_LEFT;
    char right_string[] = MATHS_CONST_RIGHT;
    long left_string_length;
    long right_string_length;
    char *result;

    left_string_length = strlen(left_string);
    right_string_length = strlen(right_string);

    while ( i < len )
    {
        switch ( (*express)[i] )
        {
            case '[':
            {
                /*
                   Numbers are kept inside square brackets, so skip this
                   as it is already (we assume) syntactically correct.
                */

                i++;

                j = len_to_sq_close((*express)+i,len-i);

                if ( j < 0 ) { return j; }

                i += j;

                break;
            }

            case 'a': case 'b': case 'c': case 'd': case 'e': case 'f':
            case 'g': case 'h': case 'i': case 'j': case 'k': case 'l':
            case 'm': case 'n': case 'o': case 'p': case 'q': case 'r':
            case 's': case 't': case 'u': case 'v': case 'w': case 'x':
            case 'y': case 'z':
            case 'A': case 'B': case 'C': case 'D': case 'E': case 'F':
            case 'G': case 'H': case 'I': case 'J': case 'K': case 'L':
            case 'M': case 'N': case 'O': case 'P': case 'Q': case 'R':
            case 'S': case 'T': case 'U': case 'V': case 'W': case 'X':
            case 'Y': case 'Z':
            case '_':
            case '*':
            case '/':
            case '+':
            case '(':
            case ')':
            case ',':
            case '^':
            {
                break;
            }

            case '0': case '1': case '2': case '3': case '4': case '5':
            case '6': case '7': case '8': case '9': case '-': case '.':
            {
                /*
                   If this is - then we need to check to ensure that this
                   is a number and not a function with unary negation
                   applied to it.  So we set j if this is a number.
                */

                j = 0;

                if ( ( (*express)[i] == '-' ) && ( j < len-1 ) )
                {
                    /*
                       Syntactically, if the - is at the end then there
                       are syntax problems here.
                    */

                    if ( i >= len-1 ) { return MERR_CONCLUDING_MINUS; }

                    switch ( (*express)[i+1] )
                    {
                        case '0':
                        case '1':
                        case '2':
                        case '3':
                        case '4':
                        case '5':
                        case '6':
                        case '7':
                        case '8':
                        case '9':
                        case '.':
                        {
                            j = 1;

                            break;
                        }

                        default:
                        {
                            break;
                        }
                    }
                }

                else
                {
                    j = 1;
                }

                if ( j )
                {
                    /*
                       First, find the end of the number.  Then enclose it
                       in appropriate brackets.  Usually, that is one[
                       on the left and ]() on the right.
                    */

                    j = len_to_real_end((*express)+i,len-i);

                    if ( j < 0 ) { return j; }

                    result = (char *) malloc((len+left_string_length+right_string_length+1)*sizeof(char));

                    if ( result == NULL )
                    {
                        return MERR_MEM_ERROR_C;
                    }

                    k = 0;

                    while ( k < i )
                    {
                        result[k] = (*express)[k];

                        k++;
                    }

                    k = 0;

                    while ( k < left_string_length )
                    {
                        result[i+k] = left_string[k];

                        k++;
                    }

                    k = 0;

                    while ( k < j )
                    {
                        result[i+left_string_length+k] = (*express)[i+k];

                        k++;
                    }

                    k = 0;

                    while ( k < right_string_length )
                    {
                        result[i+left_string_length+j+k] = right_string[k];

                        k++;
                    }

                    k = 0;

                    while ( k < len-i-j )
                    {
                        result[i+left_string_length+j+right_string_length+k] = (*express)[i+j+k];

                        k++;
                    }

                    free(*express);

                    *express = result;

                    len += left_string_length+right_string_length;

                    i += left_string_length+right_string_length+j;

                    i--;
                }

                break;
            }

            default:
            {
                return MERR_SYNTAX_ERROR_NOTHING;

                break;
            }
        }

        i++;
    }

    (*express)[len] = '\0';

    return len;
}




/*

convert_bracks - convert "bare" brackets onto proper expessions.  A bare
                 bracket can take any of the following forms:
                 - ( at the start of the function
                 - ( preceeded by any non-alpha character not including ]})

*/

long convert_bracks(char **express, long len)
{
    long i = 0;
    long j;
    long k;
    char left_string[] = MATHS_CONST_BARE_BRACK;
    long left_string_length;
    char *result;

    left_string_length = strlen(left_string);

    while ( i < len )
    {
        switch ( (*express)[i] )
        {
            case '[':
            {
                /*
                   We pass over anything in square brackets, assuming that
                   it is a number.
                */

                i++;

                j = len_to_sq_close((*express)+i,len-i);

                if ( j < 0 ) { return j; }

                i += j;

                break;
            }

            case '(':
            {
                j = 0;

                /*
                   Set j if this is a "bare" bracket.
                */

                if ( i )
                {
                    switch ( (*express)[i-1] )
                    {
                        case 'a': case 'b': case 'c': case 'd': case 'e':
                        case 'f': case 'g': case 'h': case 'i': case 'j':
                        case 'k': case 'l': case 'm': case 'n': case 'o':
                        case 'p': case 'q': case 'r': case 's': case 't':
                        case 'u': case 'v': case 'w': case 'x': case 'y':
                        case 'z':
                        case 'A': case 'B': case 'C': case 'D': case 'E':
                        case 'F': case 'G': case 'H': case 'I': case 'J':
                        case 'K': case 'L': case 'M': case 'N': case 'O':
                        case 'P': case 'Q': case 'R': case 'S': case 'T':
                        case 'U': case 'V': case 'W': case 'X': case 'Y':
                        case 'Z':
                        case '_':
                        case ')':
                        case ']':
                        {
                            break;
                        }

                        default:
                        {
                            j = 1;

                            break;
                        }
                    }
                }

                else
                {
                    j = 1;
                }

                if ( j )
                {
                    /* found a bare bracket - replace it. */

                    j = len_to_br_close((*express)+i+1,len-i-1);

                    if ( j < 0 ) { return j; }

                    j++; /* don't forget the right bracket */

                    result = (char *) malloc((len+left_string_length+1)*sizeof(char));

                    if ( result == NULL )
                    {
                        return MERR_MEM_ERROR_D;
                    }

                    k = 0;

                    while ( k < i )
                    {
                        result[k] = (*express)[k];

                        k++;
                    }

                    k = 0;

                    while ( k < left_string_length )
                    {
                        result[i+k] = left_string[k];

                        k++;
                    }

                    k = i;

                    while ( k < len )
                    {
                        result[left_string_length+k] = (*express)[k];

                        k++;
                    }

                    free(*express);

                    *express = result;

                    len += left_string_length;

                    i += left_string_length;
                }

                break;
            }

            default:
            {
                break;
            }
        }

        i++;
    }

    (*express)[len] = '\0';
    
    return len;
}




/*

convert_standbrack - standardise brackets.  All functions should be a series
                     of alphas followed by [...](...) in that order.
                     This function assumes the order is correct already, but
                     that some of the bracket pairs may be missing.

*/

long convert_standbrack(char **express, long len)
{
    long i = 0;
    long j;
    long k;
    char *result;

    while ( i < len )
    {
        switch ( (*express)[i] )
        {
            case '[':
            {
                i++;

                j = len_to_sq_close((*express)+i,len-i);

                if ( j < 0 ) { return j; }

                i += j;

                break;
            }

            case 'a': case 'b': case 'c': case 'd': case 'e': case 'f':
            case 'g': case 'h': case 'i': case 'j': case 'k': case 'l':
            case 'm': case 'n': case 'o': case 'p': case 'q': case 'r':
            case 's': case 't': case 'u': case 'v': case 'w': case 'x':
            case 'y': case 'z':
            case 'A': case 'B': case 'C': case 'D': case 'E': case 'F':
            case 'G': case 'H': case 'I': case 'J': case 'K': case 'L':
            case 'M': case 'N': case 'O': case 'P': case 'Q': case 'R':
            case 'S': case 'T': case 'U': case 'V': case 'W': case 'X':
            case 'Y': case 'Z':
            case '_':
            {
                i++;

                j = len_to_nonalpha((*express)+i,len-i);

                if ( j < 0 ) { return j; }

                i += j;

                /* OK - found a function, now see what brackets need
                   to be added. */

                /* first [] - set j if it is needed. */

                j = 0;

                if ( i < len )
                {
                    if ( (*express)[i] != '[' )
                    {
                        j = 1;
                    }
                }

                else
                {
                    j = 1;
                }

                if ( j )
                {
                    /*
                       OK, insert a [1] between the function name and
                       whatever follows.
                    */

                    result = (char *) malloc((len+4)*sizeof(char));

                    if ( result == NULL )
                    {
                        return MERR_MEM_ERROR_E;
                    }

                    k = 0;

                    while ( k < i )
                    {
                        result[k] = (*express)[k];

                        k++;
                    }

                    result[k]   = '[';
                    result[k+1] = '1';
                    result[k+2] = ']';

                    while ( k < len )
                    {
                        result[k+3] = (*express)[k];

                        k++;
                    }

                    free(*express);

                    *express = result;

                    len += 3;

                    i += 3;
                }

                else
                {
                    i++;

                    j = len_to_sq_close((*express)+i,len-i);

                    if ( j < 0 ) { return j; }

                    if ( j == 0 )
                    {
                        /*
                           [] is replaced by [1]
                        */

                        result = (char *) malloc((len+2)*sizeof(char));

                        if ( result == NULL )
                        {
                            return MERR_MEM_ERROR_F;
                        }

                        k = 0;

                        while ( k < i )
                        {
                            result[k] = (*express)[k];

                            k++;
                        }

                        result[k] = '1';

                        while ( k < len )
                        {
                            result[k+1] = (*express)[k];

                            k++;
                        }

                        free(*express);

                        *express = result;

                        len += 1;

                        i++;
                        i++;
                    }

                    else
                    {
                        i += j;

                        i++;
                    }
                }

                /* and lastly () - set j if it is missing */

                j = 0;

                if ( i < len )
                {
                    if ( (*express)[i] != '(' )
                    {
                        j = 1;
                    }
                }

                else
                {
                    j = 1;
                }

                if ( j )
                {
                    result = (char *) malloc((len+3)*sizeof(char));

                    if ( result == NULL )
                    {
                        return MERR_MEM_ERROR_G;
                    }

                    k = 0;

                    while ( k < i )
                    {
                        result[k] = (*express)[k];

                        k++;
                    }

                    result[k]   = '(';
                    result[k+1] = ')';

                    while ( k < len )
                    {
                        result[k+2] = (*express)[k];

                        k++;
                    }

                    free(*express);

                    *express = result;

                    len += 2;

                    i += 2;
                }

                break;
            }

            default:
            {
                break;
            }
        }

        i++;
    }

    (*express)[len] = '\0';

    return len;
}




/*

convert_negations  - remove unary negations by converting putting them inside
                     multiplier brackers [].

*/

long convert_negations(char **express, long len)
{
    long i = 0;
    long j;
    long k;

    while ( i < len )
    {
        switch ( (*express)[i] )
        {
            case '[':
            {
                i++;

                j = len_to_sq_close((*express)+i,len-i);

                if ( j < 0 ) { return j; }

                i += j;

                break;
            }

            case '-':
            {
                j = len_to_open((*express)+i,len-i);

                if ( j <  0 ) { return j; }
                if ( j == 0 ) { return MERR_SYNTAX_ERROR_NEG_END; }

                /*
                   Check to see if the number in the [] is -ve.  If it is
                   then the two negations will cancel.
                */

                if ( (*express)[i+j+1] == '-' )
                {
                    /* two negatives will cancel, remove both */

                    k = i;

                    while ( k < i+j )
                    {
                        (*express)[k] = (*express)[k+1];

                        k++;
                    }

                    while ( k < len-2 )
                    {
                        (*express)[k] = (*express)[k+2];

                        k++;
                    }

                    len -= 2;

                    i--;
                }

                else
                {
                    /* remove 1 negative, add another */

                    k = i;

                    while ( k < i+j )
                    {
                        (*express)[k] = (*express)[k+1];

                        k++;
                    }

                    (*express)[i+j] = '-';

                    i--;
                }

                break;
            }

            default:
            {
                break;
            }
        }

        i++;
    }

    (*express)[len] = '\0';

    return len;
}




/*

Local global constants

*/

char pow_start_expr_1_[] = MATHS_POW_STRING_LEFT;
char pow_mid_expr_1_[]   = MATHS_POW_STRING_MID;
char pow_end_expr_1_[]   = MATHS_POW_STRING_RIGHT;

char pow_sep_expr[1]    = { '^'               };
char *pow_start_expr[1] = { pow_start_expr_1_ };
char *pow_mid_expr[1]   = { pow_mid_expr_1_   };
char *pow_end_expr[1]   = { pow_end_expr_1_   };

char muldiv_start_expr_1_[] = MATHS_MUL_STRING_LEFT;
char muldiv_mid_expr_1_[]   = MATHS_MUL_STRING_MID;
char muldiv_end_expr_1_[]   = MATHS_MUL_STRING_RIGHT;

char muldiv_start_expr_2_[] = MATHS_DIV_STRING_LEFT;
char muldiv_mid_expr_2_[]   = MATHS_DIV_STRING_MID;
char muldiv_end_expr_2_[]   = MATHS_DIV_STRING_RIGHT;

char muldiv_sep_expr[2]    = { '*'                  , '/'                  };
char *muldiv_start_expr[2] = { muldiv_start_expr_1_ , muldiv_start_expr_2_ };
char *muldiv_mid_expr[2]   = { muldiv_mid_expr_1_   , muldiv_mid_expr_2_   };
char *muldiv_end_expr[2]   = { muldiv_end_expr_1_   , muldiv_end_expr_2_   };

char add_start_expr_1_[] = MATHS_ADD_STRING_LEFT;
char add_mid_expr_1_[]   = MATHS_ADD_STRING_MID;
char add_end_expr_1_[]   = MATHS_ADD_STRING_RIGHT;

char add_sep_expr[1]    = { '+'               };
char *add_start_expr[1] = { add_start_expr_1_ };
char *add_mid_expr[1]   = { add_mid_expr_1_   };
char *add_end_expr[1]   = { add_end_expr_1_   };




/*

convert_pow        - replace a^b with pow[1]{}(a,b) from right to left

*/

long convert_pow(char **express, long len)
{
    char *result;

    result = convert_rtol(1,pow_sep_expr,pow_start_expr,pow_mid_expr,pow_end_expr,*express,&len);

    if ( result != NULL )
    {
        free(*express);

        *express = result;

        len = strlen(*express);

        (*express)[len] = '\0';
    }

    return len;
}




/*

convert_muldiv     - replace a*b with mul[1]{}(a,b) from left to right and
                     a/b with mul[1]{}(a,inv(b)) left to right simultaneously

*/

long convert_muldiv(char **express, long len)
{
    char *result;

    result = convert_ltor(2,muldiv_sep_expr,muldiv_start_expr,muldiv_mid_expr,muldiv_end_expr,*express,&len);

    if ( result != NULL )
    {
        free(*express);

        *express = result;

        len = strlen(*express);

        (*express)[len] = '\0';
    }

    return len;
}




/*

convert_add        - replace a+b with add[1]{}(a,b) from left to right

*/

long convert_add(char **express, long len)
{
    char *result;

    result = convert_ltor(1,add_sep_expr,add_start_expr,add_mid_expr,add_end_expr,*express,&len);

    if ( result != NULL )
    {
        free(*express);

        *express = result;

        len = strlen(*express);

        (*express)[len] = '\0';
    }

    return len;
}








char *convert_ltor(int lev_size, char *sep_expr, char **start_expr, char **mid_expr, char **end_expr, char *express, long *len)
{
    long i,j,k,jb,kb;
    long q;
    char *tempstr;
    char *startstr;
    char *startstrb;
    char *startstrc;
    char *result;
    char *whatb;
    long b_zero;
    long lenb;
    int condind;

    /* FIXME: this needs documenting */

    /* Find the number of characters before either the brackets open
       or the expression ends. */

    j = len_to_br_open(express,(*len));

    if ( j < 0 ) { (*len) = j; return NULL; }

    /* Find the number of characters taken by the expression (...),
       assuming that such exists. */

    k = 2 + len_to_br_close(express+j+1,(*len)-j-1);

    if ( k < 2 ) { (*len) = k-2; return NULL; }

    /* j+k is the total length of the expression (including brackets and
       express is in them).
       k-2 is the length of the expression in the brackets
      
       eg. sin[2](x) => j+k = 9  sin[2](x)    &(express[0])
                        j = 6    sin[2]       &(express[0])
                        k = 3    (x)          &(express[j=6])
                        k-2 = 1  x            &(express[j+1=7])
      
       Next character is at express[j+k] */

    /* Work out this expression. */

    if ( k > 2 )
    {
        q = k-2;

        tempstr = convert_ltor(lev_size,sep_expr,start_expr,mid_expr,end_expr,express+j+1,&q);

        if ( tempstr == NULL )
        {
            *len = q;

            return NULL;
        }

        startstr = (char *) malloc((j+1+strlen(tempstr)+2)*sizeof(char));

        if ( startstr == NULL )
        {
            free(tempstr);

            *len = MERR_MEM_ERRORB_A;

            return NULL;
        }

        strncpy(startstr,express,j+1);
        startstr[j+1] = '\0';
        strcat(startstr,tempstr);
        strcat(startstr,")");

        free(tempstr);
    }

    else
    {
        startstr = (char *) malloc((j+3)*sizeof(char));

        if ( startstr == NULL )
        {
            (*len) = MERR_MEM_ERRORB_B;

            return NULL;
        }

        strncpy(startstr,express,j+2);
        startstr[j+2] = '\0';
    }

    /* If this is the end of the expression then we just copy this and
       return. */

    if ( (*len) == j+k )
    {
        result = startstr;
    }

    else
    {
        /* If the next separator is not a separator, then recurse. */

        condind = 0;

        for ( i = 1 ; i <= lev_size ; i++ )
        {
            if ( express[j+k] == sep_expr[i-1] )
            {
                condind = i;
            }
        }

        if ( condind )
        {
            b_zero = j+k+1;
            lenb   = (*len) - b_zero;
            whatb  = express+b_zero;

            jb = len_to_br_open(whatb,lenb);

            if ( jb < 0 ) { (*len) = jb; free(startstr); return NULL; }

            kb = 2 + len_to_br_close(whatb+jb+1,lenb-jb-1);

            if ( kb < 2 ) { (*len) = kb-2; free(startstr); return NULL; }

            if ( kb > 2 )
            {
                q = kb-2;

                tempstr = convert_ltor(lev_size,sep_expr,start_expr,mid_expr,end_expr,whatb+jb+1,&q);

                if ( tempstr == NULL )
                {
                    *len = q;

                    free(startstr);

                    return NULL;
                }

                startstrb = (char *) malloc((jb+1+strlen(tempstr)+2)*sizeof(char));

                if ( startstrb == NULL )
                {
                    free(tempstr);
                    free(startstr);

                    (*len) = MERR_MEM_ERRORB_C;

                    return NULL;
                }

                strncpy(startstrb,whatb,jb+1);
                startstrb[jb+1] = '\0';
                strcat(startstrb,tempstr);
                strcat(startstrb,")");

                free(tempstr);
            }

            else
            {
                startstrb = (char *) malloc((jb+3)*sizeof(char));

                if ( startstrb == NULL )
                {
                    free(startstr);

                    (*len) = MERR_MEM_ERRORB_D;

                    return NULL;
                }

                strncpy(startstrb,whatb,jb+2);
                startstrb[jb+2] = '\0';
            }

            /* Now we make the expression */

            startstrc = (char *) malloc((strlen(start_expr[condind-1])+strlen(startstr)+strlen(mid_expr[condind-1])+strlen(startstrb)+strlen(end_expr[condind-1])+1)*sizeof(char));

            if ( startstrc == NULL )
            {
                free(startstr);
                free(startstrb);

                (*len) = MERR_MEM_ERRORB_E;

                return NULL;
            }

            strcpy(startstrc,start_expr[condind-1]);
            strcat(startstrc,startstr);
            strcat(startstrc,mid_expr[condind-1]);
            strcat(startstrc,startstrb);
            strcat(startstrc,end_expr[condind-1]);

            free(startstr);
            free(startstrb);

            /* And put it all together */

            result = (char *) malloc((strlen(startstrc)+(lenb-(jb+kb))+1)*sizeof(char));

            if ( result == NULL )
            {
                free(startstrc);

                (*len) = MERR_MEM_ERRORB_F;

                return NULL;
            }

            strcpy(result,startstrc);
            strncat(result,whatb+jb+kb,lenb-jb-kb);
            result[strlen(startstrc)+(lenb-(jb+kb))] = '\0';

            free(startstrc);

            /* do others */

            startstrc = result;

            q = strlen(startstrc);

            result = convert_ltor(lev_size,sep_expr,start_expr,mid_expr,end_expr,startstrc,&q);

            if ( result == NULL )
            {
                *len = q;

                return NULL;
            }

            free(startstrc);
        }

        else
        {
            q = (*len)-j-k-1;

            tempstr = convert_ltor(lev_size,sep_expr,start_expr,mid_expr,end_expr,express+j+k+1,&q);

            if ( tempstr == NULL )
            {
                *len = q;

                free(startstr);

                return NULL;
            }

            result = (char *) malloc((strlen(startstr)+1+strlen(tempstr)+1)*sizeof(char));

            if ( result == NULL )
            {
                free(tempstr);
                free(startstr);

                *len = MERR_MEM_ERRORB_M;

                return NULL;
            }

            strcpy(result,startstr);
            strncat(result,express+j+k,1);
            result[strlen(startstr)+1] = '\0';
            strcat(result,tempstr);

            free(tempstr);
            free(startstr);
        }
    }

    (*len) = strlen(result);

    return result;
}

char *convert_rtol(int lev_size, char *sep_expr, char **start_expr, char **mid_expr, char **end_expr, char *express, long *len)
{
    long i,j,k,jb,kb;
    long q;
    char *tempstr;
    char *startstr;
    char *startstrb;
    char *startstrc;
    char *result;
    char *endstr;
    long lenb;
    int condind;

    /* FIXME: this needs documenting */

    /* Find the number of characters before either the brackets open
       or the expression ends. */

    j = len_to_br_open(express,(*len));

    if ( j < 0 ) { (*len) = j; return NULL; }

    /* Find the number of characters taken by the expression (...),
       assuming that such exists. */

    k = 2 + len_to_br_close(express+j+1,(*len)-j-1);

    if ( k < 2 ) { (*len) = k-2; return NULL; }

    /* j+k is the total length of the expression (including brackets and
       express is in them).
       k-2 is the length of the expression in the brackets
      
       eg. sin[2](x) => j+k = 9  sin[2](x)    &(express[0])
                        j = 6    sin[2]       &(express[0])
                        k = 3    (x)          &(express[j=6])
                        k-2 = 1  x            &(express[j+1=7])
      
       Next character is at express[j+k] */

    /* Work out this expression */

    if ( k > 2 )
    {
        q = k-2;

        tempstr = convert_rtol(lev_size,sep_expr,start_expr,mid_expr,end_expr,express+j+1,&q);

        if ( tempstr == NULL )
        {
            *len = q;

            return NULL;
        }

        startstr = (char *) malloc((j+1+strlen(tempstr)+2)*sizeof(char));

        if ( startstr == NULL )
        {
            free(tempstr);

            (*len) = MERR_MEM_ERRORB_G;

            return NULL;
        }

        strncpy(startstr,express,j+1);
        startstr[j+1] = '\0';
        strcat(startstr,tempstr);
        strcat(startstr,")");

        free(tempstr);
    }

    else
    {
        startstr = (char *) malloc((j+3)*sizeof(char));

        if ( startstr == NULL )
        {
            (*len) = MERR_MEM_ERRORB_H;

            return NULL;
        }

        strncpy(startstr,express,j+2);
        startstr[j+2] = '\0';
    }

    /* If this is the end of the expression then we just copy this and
       return. */

    if ( (*len) == j+k )
    {
        result = startstr;
    }

    else
    {
        /* Recurse the rest of the string before checking here. */

        q = (*len)-j-k-1;

        endstr = convert_rtol(lev_size,sep_expr,start_expr,mid_expr,end_expr,express+j+k+1,&q);

        if ( endstr == NULL )
        {
            *len = q;

            free(startstr);

            return NULL;
        }

        /* So we have:
           "startstr"."_"."endstr"
           where _ is some operator */

        /* If we have a multiply, then fix it. */

        condind = 0;

        for ( i = 1 ; i <= lev_size ; i++ )
        {
            if ( express[j+k] == sep_expr[i-1] )
            {
                condind = i;
            }
        }

        if ( condind )
        {
            /* Figure out expresss to the right. */

            lenb = strlen(endstr);

            jb = len_to_br_open(endstr,lenb);

            if ( jb < 0 ) { (*len) = jb; free(startstr); free(endstr); return NULL; }

            kb = 2 + len_to_br_close(endstr+jb+1,lenb-jb-1);

            if ( kb < 2 ) { (*len) = kb-2; free(startstr); free(endstr); return NULL; }

            /* Don't worry about express's in brackets on the right - it will
               have been done already anyhow. */

            startstrb = (char *) malloc((jb+kb+1)*sizeof(char));

            if ( startstrb == NULL )
            {
                free(startstr);
                free(endstr);

                (*len) = MERR_MEM_ERRORB_I;

                return NULL;
            }

            strncpy(startstrb,endstr,jb+kb);
            startstrb[jb+kb] = '\0';

            /* Now we make the expression */

            startstrc = (char *) malloc((strlen(start_expr[condind-1])+strlen(startstr)+strlen(mid_expr[condind-1])+strlen(startstrb)+strlen(end_expr[condind-1])+1)*sizeof(char));

            if ( startstrc == NULL )
            {
                free(startstr);
                free(startstrb);
                free(endstr);

                (*len) = MERR_MEM_ERRORB_J;

                return NULL;
            }

            strcpy(startstrc,start_expr[condind-1]);
            strcat(startstrc,startstr);
            strcat(startstrc,mid_expr[condind-1]);
            strcat(startstrc,startstrb);
            strcat(startstrc,end_expr[condind-1]);

            free(startstr);
            free(startstrb);

            /* And put it all together */

            result = (char *) malloc((strlen(startstrc)+(lenb-(jb+kb))+1)*sizeof(char));

            if ( result == NULL )
            {
                free(startstrc);
                free(endstr);

                (*len) = MERR_MEM_ERRORB_K;

                return NULL;
            }

            strcpy(result,startstrc);
            strncat(result,endstr+jb+kb,lenb-jb-kb);
            result[strlen(startstrc)+(lenb-(jb+kb))] = '\0';

            free(startstrc);
            free(endstr);
        }

        else
        {
            result = (char *) malloc((strlen(startstr)+1+strlen(endstr)+1)*sizeof(char));

            if ( result == NULL )
            {
                free(startstr);
                free(endstr);

                (*len) = MERR_MEM_ERRORB_L;

                return NULL;
            }

            strcpy(result,startstr);
            result[strlen(startstr)] = express[j+k];
            result[strlen(startstr)+1] = '\0';
            strcat(result,endstr);

            free(startstr);
            free(endstr);
        }
    }

    (*len) = strlen(result);

    return result;
}




/*

convert_passes     - replace expressions pass[1](...) with ...

*/

long convert_passes(char **express, long len)
{
    long i = 0;
    long j,k,l;

    while ( i < len )
    {
        switch ( (*express)[i] )
        {
            case '(':
            case ')':
            case ',':
            {
                break;
            }

            case '[':
            {
                i++;

                j = len_to_sq_close((*express)+i,len-i);

                if ( j <  0 ) { return j; }
                if ( j == 0 ) { return MERR_SYNTAX_ERROR_PASSREM_A; }

                i += j;

                break;
            }

            case 'a': case 'b': case 'c': case 'd': case 'e': case 'f':
            case 'g': case 'h': case 'i': case 'j': case 'k': case 'l':
            case 'm': case 'n': case 'o': case 'p': case 'q': case 'r':
            case 's': case 't': case 'u': case 'v': case 'w': case 'x':
            case 'y': case 'z':
            case 'A': case 'B': case 'C': case 'D': case 'E': case 'F':
            case 'G': case 'H': case 'I': case 'J': case 'K': case 'L':
            case 'M': case 'N': case 'O': case 'P': case 'Q': case 'R':
            case 'S': case 'T': case 'U': case 'V': case 'W': case 'X':
            case 'Y': case 'Z':
            case '_':
            {
                j = len_to_open((*express)+i,len-i);

                if ( j < 0 ) { return j; }

                if ( i+j+1 >= len ) { return MERR_SYNTAX_ERROR_PASSREM_O; }

                if ( ( strlen(MATHS_CONST_BARE_BRACK) == j                 ) &&
                     ( strncmp((*express)+i,MATHS_CONST_BARE_BRACK,j) == 0 )    )
                {
                    switch ( (*express)[i+j+1] )
                    {
                        case '1':
                        {
                            i += j;

                            if ( i >= len )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_O;
                            }

                            if ( (*express)[i] != '[' )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_C;
                            }

                            i++;

                            if ( i >= len )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_O;
                            }

                            if ( (*express)[i] != '1' )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_D;
                            }

                            i++;

                            if ( i >= len )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_O;
                            }

                            if ( (*express)[i] != ']' )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_E;
                            }

                            i++;

                            if ( i >= len )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_O;
                            }

                            if ( (*express)[i] != '(' )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_F;
                            }

                            i++;

                            break;
                        }

                        case '-':
                        {
                            i += j;

                            if ( i >= len )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_O;
                            }

                            if ( (*express)[i] != '[' )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_G;
                            }

                            i++;

                            if ( i >= len )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_O;
                            }

                            if ( (*express)[i] != '-' )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_H;
                            }

                            i++;

                            if ( i >= len )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_O;
                            }

                            if ( (*express)[i] != '1' )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_I;
                            }

                            i++;

                            if ( i >= len )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_O;
                            }

                            if ( (*express)[i] != ']' )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_J;
                            }

                            i++;

                            if ( i >= len )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_O;
                            }

                            if ( (*express)[i] != '(' )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_K;
                            }

                            i++;

                            k = len_to_open((*express)+i,len-i);

                            if ( i+k >= len )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_O;
                            }

                            if ( (*express)[i+k] != '[' )
                            {
                                return MERR_SYNTAX_ERROR_PASSREM_N;
                            }

                            l = i-4;

                            while ( l < i+k )
                            {
                                (*express)[l] = (*express)[l+1];

                                l++;
                            }

                            (*express)[l] = '-';

                            i--;

                            break;
                        }

                        default:
                        {
                            return MERR_SYNTAX_ERROR_PASSREM_L;

                            break;
                        }
                    }

                    k = len_to_br_close((*express)+i,len-i);

                    j += 4;
                    i -= j;

                    l = i+j+k+1;

                    while ( l < len )
                    {
                        (*express)[l-1] = (*express)[l];

                        l++;
                    }

                    len--;

                    l = i+j;

                    while ( l < len )
                    {
                        (*express)[l-j] = (*express)[l];

                        l++;
                    }

                    len -= j;
                }

                else
                {
                    i += j;

                    i--;
                }

                break;
            }

            default:
            {
                return MERR_SYNTAX_ERROR_PASSREM_M;

                break;
            }
        }

        i++;
    }

    (*express)[len] = '\0';

    return len;
}

