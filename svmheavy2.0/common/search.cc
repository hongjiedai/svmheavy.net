
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
// Stream utility
//
// Written by: Alistair Shilton
//             Melbourne University
//

//
// Stream clearing - runs through until a ':' is encountered
//


#include <iostream>
#include "search.h"


//
// off-load stream until after next ': '
//

std::istream &operator>>(std::istream &input, wait_dummy &)
{
    char scanner = 'z';

    while ( scanner != ':' )
    {
        input >> scanner;
    }

    input.ignore(1);

    return input;
}



