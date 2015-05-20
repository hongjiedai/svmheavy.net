
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
// Support vector machine optimisation functions
//
// Written by: Alistair Shilton
//             Melbourne University
//

#ifndef _svoptim_h
#define _svoptim_h

#include "svflags.h"
#include "svdata.h"
#include "kernel.h"
#include "common/vector.h"
#include "common/matrix.h"
#include "common/factor.h"

#ifdef DO__FLOPS
#define JTYPE L_DOUBLE
#endif
#ifndef DO__FLOPS
#define JTYPE volatile L_DOUBLE
#endif

unsigned long solve(SVdata &problem, L_DOUBLE &sol_tol, unsigned long epochs, volatile int *async_exit_flag, volatile unsigned long *loop_count, volatile unsigned long *loop_count_last, JTYPE *J);

#endif
