
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
// Thread abstraction
//
// Written by: Alistair Shilton
//             Melbourne University
//


#ifndef _svm_thread_h
#define _svm_thread_h

/*
   svmthread_create - creates a detached thread which runs function
                      start_routine(arg).  Will return 0 on success, or
                      nonzero (error code) on failure.  In djgpp will
                      always return -1 (failure), as dos cannot do
                      threads.

   svmthread_test - return nonzero if threads available
*/

int svmthread_create(void * (*start_routine)(void *), void *arg);
int svmthread_test(void);

#endif
