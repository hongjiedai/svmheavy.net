#include "stdafx.h"
/*
 *  RIMElib: RuntIme Mathematical Equation Library (example code)
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


#define TOP_LEVEL

#include "rime.h"
#include <stdio.h>
#include <iostream>


double x_global;


double DoubleArgEvaluationFunction(void *argContents, long i, long j);


int main()
{
    char sprock[512];
    MA_node *what;
    MA_node *whet;
    int physics_colloqium;
    long count = 0;
    long n,m,i,j,k,l;
    double x,y,z;

    repeat:

    std::cout << "Equation: "; std::cin >> sprock;

    try { what = new MA_node(sprock); }

    catch ( int x )
    {
        std::cout << "Error parsing equation " << x << "\n";

        exit(1);
    }

    again:

    count++;

    std::cout << "count = " << count << "\n";

    std::cout << "1. Print expression.\n";
    std::cout << "2. Differentiate.\n";
    std::cout << "3. Reenter expression.\n";
    std::cout << "4. Evaluate for given x.\n";
    std::cout << "5. Evaluate for given x,y.\n";
    std::cout << "6. Evaluate for given x,y,z.\n";
    std::cout << "7. Test input stream.\n";
    std::cout << "8. Substitute variable.\n";
    std::cout << "9. Row substitute variable.\n";
    std::cout << "10. Column substitute variable.\n";
    std::cout << "11. Test general evaluation.\n";
    std::cout << "12. Matrix Differentiate.\n";
    std::cout << "13. Row Differentiate.\n";
    std::cout << "14. Column Differentiate.\n";

    std::cin >> physics_colloqium;

    switch ( physics_colloqium )
    {
        case 1:
        {
            std::cout << *what;

            std::cout << "\n";

            break;
        }

        case 2:
        {
            std::cout << "n: "; std::cin >> n;
            std::cout << "m: "; std::cin >> m;

            whet = what;

            try { what = what->diffp(n,m); }

            catch ( int x )
            {
                std::cout << "Error differentiating " << x << "\n";

                exit(1);
            }

            delete whet;

            break;
        }

        case 3:
        {
            goto repeat;

            break;
        }

        case 4:
        {
            std::cout << "x = "; std::cin >> x;

            std::cout << "result = " << (*what)(x) << "\n";

            break;
        }

        case 5:
        {
            std::cout << "x = "; std::cin >> x;
            std::cout << "y = "; std::cin >> y;

            std::cout << "result = " << (*what)(x,y) << "\n";

            break;
        }

        case 6:
        {
            std::cout << "x = "; std::cin >> x;
            std::cout << "y = "; std::cin >> y;
            std::cout << "z = "; std::cin >> z;

            std::cout << "result = " << (*what)(x,y,z) << "\n";

            break;
        }

        case 7:
        {
            delete what;

            what = new MA_node;

            try
            {
                std::cin >> *what;
            }

            catch ( int x )
            {
                printf("Error %d during stream operation.\n",x);

                return 1;
            }

            break;
        }

        case 8:
        {
            std::cout << "n = "; std::cin >> n;
            std::cout << "m = "; std::cin >> m;

            whet = new MA_node;

            std::cout << "var(" << n << "," << m << ") = "; std::cin >> *whet;

            what->replace(*whet,n,m);

            delete whet;

            break;
        }

        case 9:
        {
            std::cout << "n = "; std::cin >> n;

            whet = new MA_node;

            std::cout << "var(" << n << ",var(0,0)) = "; std::cin >> *whet;

            what->replaceRow(*whet,n);

            delete whet;

            break;
        }

        case 10:
        {
            std::cout << "m = "; std::cin >> m;

            whet = new MA_node;

            std::cout << "var(var(0,0)," << m << ") = "; std::cin >> *whet;

            what->replaceCol(*whet,m);

            delete whet;

            break;
        }

        case 11:
        {
            std::cout << "x = "; std::cin >> x_global;

            std::cout << "result = " << (*what)(DoubleArgEvaluationFunction,&x_global) << "\n";

            break;
        }

        case 12:
        {
            std::cout << "i: "; std::cin >> i;
            std::cout << "j: "; std::cin >> j;
            std::cout << "k: "; std::cin >> k;
            std::cout << "l: "; std::cin >> l;

            whet = what;

            try { what = what->diffpMatrix(i,j,k,l); }

            catch ( int x )
            {
                std::cout << "Error differentiating " << x << "\n";

                exit(1);
            }

            delete whet;

            break;
        }

        case 13:
        {
            std::cout << "n: "; std::cin >> n;
            std::cout << "k: "; std::cin >> k;
            std::cout << "l: "; std::cin >> l;

            whet = what;

            try { what = what->diffpRow(n,k,l); }

            catch ( int x )
            {
                std::cout << "Error differentiating " << x << "\n";

                exit(1);
            }

            delete whet;

            break;
        }

        case 14:
        {
            std::cout << "i: "; std::cin >> i;
            std::cout << "j: "; std::cin >> j;
            std::cout << "m: "; std::cin >> m;

            whet = what;

            try { what = what->diffpCol(i,j,m); }

            catch ( int x )
            {
                std::cout << "Error differentiating " << x << "\n";

                exit(1);
            }

            delete whet;

            break;
        }

        default:
        {
            delete what;

            return 0;

            break;
        }
    }

    goto again;

    return 0;
}




double DoubleArgEvaluationFunction(void *argContents, long i, long j)
{
    std::cout << "var(" << i << "," << j << ") = " << *((double *) argContents) << "\n";

    return *((double *) argContents);
}
