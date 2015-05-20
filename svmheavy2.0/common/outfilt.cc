
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


//
// Output filter (removes inf and nan)
//
// Written by: Alistair Shilton
//             Melbourne University
//

#include <float.h>
#include "c_double.h"

char *fixnum(double input)
{
    static char buffer[128];

    sprintf(buffer,"%.12le",input);

    if ( !( ( buffer[0] == '0' ) ||
            ( buffer[0] == '1' ) ||
            ( buffer[0] == '2' ) ||
            ( buffer[0] == '3' ) ||
            ( buffer[0] == '4' ) ||
            ( buffer[0] == '5' ) ||
            ( buffer[0] == '6' ) ||
            ( buffer[0] == '7' ) ||
            ( buffer[0] == '8' ) ||
            ( buffer[0] == '9' ) ||
            ( buffer[0] == '.' ) ||
            ( buffer[0] == '-' ) ||
            ( buffer[0] == '+' )    ) )
    {
        if ( input < 0 )
        {
            if ( input >= 0 )
            {
                sprintf(buffer,"              0.000");
            }

            else
            {
                sprintf(buffer,"             %.12le",DBL_MAX);
            }
        }

        else if ( input > 0 )
        {
            if ( input <= 0 )
            {
                sprintf(buffer,"              0.000");
            }

            else
            {
                sprintf(buffer,"             %.12le",DBL_MIN);
            }
        }

        else
        {
            sprintf(buffer,"              0.000");
        }
    }

    return buffer;
}

char *fixnum(c_double input)
{
    double temp = 0;

    temp += input;

    return fixnum(temp);
}


