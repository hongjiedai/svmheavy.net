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
#include "stdafx.h"
/* 先註解起來以便編譯成功

//
// Noise generation example.
//


#define TOP_LEVEL

#include "rime.h"
#include <stdio.h>
#include <iostream>


double x_global;


int main()
{
    long i,n,d;
    double x,ave,stdev;

    MA_node genfn("prand(0,1,x)");

    std::cerr << "Number of points: "; std::cin >> n;
    std::cerr << "Degree: "; std::cin >> d;

    ave = 0;
    stdev = 0;

    for ( i = 0 ; i < n ; i++ )
    {
        x = genfn(d);

	ave += x/((double) n);
	stdev += (x*x)/((double) n);

	std::cout << x << "\n";
    }

    std::cerr << "Average = " << ave << "\n";
    std::cerr << "Variance = " << stdev << "\n";

    return 0;
}
*/