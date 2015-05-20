
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
// Levi civita function.
//
// Written by: Alistair Shilton
//             Melbourne University
//

#ifndef _vecmatbase_h
#define _vecmatbase_h

//
// levi_civita: calculate levi-civita terms.
//              dim = # terms in expansion (>= 1)
//              terms = some permutation of {1,2,...,dim}
//                      (or a subset repeated to make dim terms)
//

int levi_civita(long dim, long *terms);

#endif
