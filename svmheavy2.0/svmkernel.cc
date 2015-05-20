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
// Support vector machine - kernel file creation utility
//
// Written by: Alistair Shilton
//             Melbourne University
//


#define TOP_LEVEL

#include <iostream>
#include <fstream>
#include "nullstream.h"
#include "gsl/err/gsl_errno.h"
#include <time.h>
#include "common/svdefs.h"
#include "common/sparsevector.h"
#include "svm_pattern.h"
#include "svm_regress.h"


MEM_DEFS()

int main()
{
    try
    {
        std::cout << "SVMkernel 1.1 - an SVM kernel creator for SVMheavy by Alistair Shilton.\n";
        std::cout << "=============\n\n";

        char buffer[1024];

        std::cout << "Kernel filename: "; std::cin >> buffer;
        std::ofstream kerfile(buffer);

        if ( !(kerfile.is_open()) )
        {
            std::cerr << "Unable to open file\n";
            exit(1);
        }

        kernel_dicto K;

        std::cout << "Enter kernel details (sparsity must match training data).\n\n";

        K.prime_io(&(std::cout),0);

        std::cin >> K;

        kerfile << K;

        kerfile.close();
    }

    catch ( int x )
    {
        REPORT_THROW(std::cout,x);
    }

    MEM_REPORT();

    return 0;
}
