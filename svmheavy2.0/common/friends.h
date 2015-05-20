
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


#ifndef _friends_h
#define _friends_h

#define ALL_FRIENDS                                                                     \
friend std::ostream &operator<<(std::ostream &output, const fVECTOR   &source);         \
friend std::ostream &operator<<(std::ostream &output, const fMATRIX   &source);         \
friend std::ostream &operator<<(std::ostream &output, const fMATRFact &source);         \
                                                                                        \
friend std::istream &operator>>(std::istream &input, fVECTOR   &dest);                  \
friend std::istream &operator>>(std::istream &input, fMATRIX   &dest);                  \
friend std::istream &operator>>(std::istream &input, fMATRFact &dest);                  \
                                                                                        \
friend fVECTOR   operator+ (const fVECTOR  &left_op);                                   \
friend fMATRIX   operator+ (const fMATRIX  &left_op);                                   \
                                                                                        \
friend fVECTOR   operator- (const fVECTOR  &left_op);                                   \
friend fMATRIX   operator- (const fMATRIX  &left_op);                                   \
                                                                                        \
friend fVECTOR   operator+ (const fVECTOR  &left_op, const fVECTOR  &right_op);         \
friend fVECTOR   operator+ (const fVECTOR  &left_op, const   double &right_op);         \
friend fVECTOR   operator+ (const fVECTOR  &left_op, const c_double &right_op);         \
friend fVECTOR   operator+ (const   double &left_op, const fVECTOR  &right_op);         \
friend fVECTOR   operator+ (const c_double &left_op, const fVECTOR  &right_op);         \
friend fMATRIX   operator+ (const fMATRIX  &left_op, const fMATRIX  &right_op);         \
                                                                                        \
friend fVECTOR   operator- (const fVECTOR  &left_op, const fVECTOR  &right_op);         \
friend fVECTOR   operator- (const fVECTOR  &left_op, const   double &right_op);         \
friend fVECTOR   operator- (const fVECTOR  &left_op, const c_double &right_op);         \
friend fVECTOR   operator- (const   double &left_op, const fVECTOR  &right_op);         \
friend fVECTOR   operator- (const c_double &left_op, const fVECTOR  &right_op);         \
friend fMATRIX   operator- (const fMATRIX  &left_op, const fMATRIX  &right_op);         \
                                                                                        \
friend fVECTOR  &operator+=(      fVECTOR  &left_op, const fVECTOR  &right_op);         \
friend fVECTOR  &operator+=(      fVECTOR  &left_op, const   double &right_op);         \
friend fVECTOR  &operator+=(      fVECTOR  &left_op, const c_double &right_op);         \
friend fMATRIX  &operator+=(      fMATRIX  &left_op, const fMATRIX  &right_op);         \
                                                                                        \
friend fVECTOR  &operator-=(      fVECTOR  &left_op, const fVECTOR  &right_op);         \
friend fVECTOR  &operator-=(      fVECTOR  &left_op, const   double &right_op);         \
friend fVECTOR  &operator-=(      fVECTOR  &left_op, const c_double &right_op);         \
friend fMATRIX  &operator-=(      fMATRIX  &left_op, const fMATRIX  &right_op);         \
                                                                                        \
friend fVECTOR   operator* (const fVECTOR  &left_op, const   double &right_op);         \
friend fVECTOR   operator* (const fVECTOR  &left_op, const c_double &right_op);         \
friend fMATRIX   operator* (const fMATRIX  &left_op, const   double &right_op);         \
friend fMATRIX   operator* (const fMATRIX  &left_op, const c_double &right_op);         \
                                                                                        \
friend fVECTOR   operator* (const   double &left_op, const fVECTOR  &right_op);         \
friend fVECTOR   operator* (const c_double &left_op, const fVECTOR  &right_op);         \
friend fMATRIX   operator* (const   double &left_op, const fMATRIX  &right_op);         \
friend fMATRIX   operator* (const c_double &left_op, const fMATRIX  &right_op);         \
                                                                                        \
friend fVECTOR   operator/ (const fVECTOR  &left_op, const   double &right_op);         \
friend fVECTOR   operator/ (const fVECTOR  &left_op, const c_double &right_op);         \
friend fMATRIX   operator/ (const fMATRIX  &left_op, const   double &right_op);         \
friend fMATRIX   operator/ (const fMATRIX  &left_op, const c_double &right_op);         \
                                                                                        \
friend L_DOUBLE  operator* (const fVECTOR  &left_op, const fVECTOR  &right_op);         \
friend fMATRIX   operator* (const fMATRIX  &left_op, const fMATRIX  &right_op);         \
friend fVECTOR   operator* (const fMATRIX  &left_op, const fVECTOR  &right_op);         \
friend fVECTOR   operator* (const fVECTOR  &left_op, const fMATRIX  &right_op);         \
                                                                                        \
friend fVECTOR  &operator*=(      fVECTOR  &left_op, const   double &right_op);         \
friend fVECTOR  &operator*=(      fVECTOR  &left_op, const c_double &right_op);         \
friend fMATRIX  &operator*=(      fMATRIX  &left_op, const   double &right_op);         \
friend fMATRIX  &operator*=(      fMATRIX  &left_op, const c_double &right_op);         \
                                                                                        \
friend fMATRIX  &assign_all(      fMATRIX  &left_op, const   double &right_op);         \
friend fMATRIX  &assign_all(      fMATRIX  &left_op, const c_double &right_op);         \
                                                                                        \
friend fVECTOR  &operator/=(      fVECTOR  &left_op, const   double &right_op);         \
friend fVECTOR  &operator/=(      fVECTOR  &left_op, const c_double &right_op);         \
friend fMATRIX  &operator/=(      fMATRIX  &left_op, const   double &right_op);         \
friend fMATRIX  &operator/=(      fMATRIX  &left_op, const c_double &right_op);         \
                                                                                        \
friend fMATRIX  &operator*=(      fMATRIX  &left_op, const fMATRIX  &right_op);         \
friend fVECTOR  &operator*=(      fVECTOR  &left_op, const fMATRIX  &right_op);         \
                                                                                        \
friend int operator==(const fVECTOR  &left_op, const fVECTOR  &right_op);               \
friend int operator==(const fMATRIX  &left_op, const fMATRIX  &right_op);               \
                                                                                        \
friend int operator==(const fVECTOR  &left_op, const   double &right_op);               \
friend int operator==(const fVECTOR  &left_op, const c_double &right_op);               \
friend int operator==(const   double &left_op, const fVECTOR  &right_op);               \
friend int operator==(const c_double &left_op, const fVECTOR  &right_op);               \
                                                                                        \
friend int operator!=(const fVECTOR  &left_op, const fVECTOR  &right_op);               \
friend int operator!=(const fMATRIX  &left_op, const fMATRIX  &right_op);               \
                                                                                        \
friend int operator!=(const fVECTOR  &left_op, const   double &right_op);               \
friend int operator!=(const fVECTOR  &left_op, const c_double &right_op);               \
friend int operator!=(const   double &left_op, const fVECTOR  &right_op);               \
friend int operator!=(const c_double &left_op, const fVECTOR  &right_op);               \
                                                                                        \
friend int operator< (const fVECTOR  &left_op, const fVECTOR  &right_op);               \
friend int operator< (const fMATRIX  &left_op, const fMATRIX  &right_op);               \
                                                                                        \
friend int operator< (const fVECTOR  &left_op, const   double &right_op);               \
friend int operator< (const fVECTOR  &left_op, const c_double &right_op);               \
friend int operator< (const   double &left_op, const fVECTOR  &right_op);               \
friend int operator< (const c_double &left_op, const fVECTOR  &right_op);               \
                                                                                        \
friend int operator<=(const fVECTOR  &left_op, const fVECTOR  &right_op);               \
friend int operator<=(const fMATRIX  &left_op, const fMATRIX  &right_op);               \
                                                                                        \
friend int operator<=(const fVECTOR  &left_op, const   double &right_op);               \
friend int operator<=(const fVECTOR  &left_op, const c_double &right_op);               \
friend int operator<=(const   double &left_op, const fVECTOR  &right_op);               \
friend int operator<=(const c_double &left_op, const fVECTOR  &right_op);               \
                                                                                        \
friend int operator> (const fVECTOR  &left_op, const fVECTOR  &right_op);               \
friend int operator> (const fMATRIX  &left_op, const fMATRIX  &right_op);               \
                                                                                        \
friend int operator> (const fVECTOR  &left_op, const   double &right_op);               \
friend int operator> (const fVECTOR  &left_op, const c_double &right_op);               \
friend int operator> (const   double &left_op, const fVECTOR  &right_op);               \
friend int operator> (const c_double &left_op, const fVECTOR  &right_op);               \
                                                                                        \
friend int operator>=(const fVECTOR  &left_op, const fVECTOR  &right_op);               \
friend int operator>=(const fMATRIX  &left_op, const fMATRIX  &right_op);               \
                                                                                        \
friend int operator>=(const fVECTOR  &left_op, const   double &right_op);               \
friend int operator>=(const fVECTOR  &left_op, const c_double &right_op);               \
friend int operator>=(const   double &left_op, const fVECTOR  &right_op);               \
friend int operator>=(const c_double &left_op, const fVECTOR  &right_op);               \
                                                                                        \
friend std::ostream &operator<<(std::ostream &output, const slow_f_fVECTOR &source);    \
friend std::istream &operator>>(std::istream &input,        slow_f_fVECTOR &dest  );    \
friend std::ostream &operator<<(std::ostream &output, const f_fVECTOR &source);         \
friend std::istream &operator>>(std::istream &input,        f_fVECTOR &dest  );         \
friend std::ostream &operator<<(std::ostream &output, const vpVECTOR &source);          \
friend std::istream &operator>>(std::istream &input,        vpVECTOR &dest  );          \
friend std::ostream &operator<<(std::ostream &output, const iVECTOR   &source);         \
friend std::istream &operator>>(std::istream &input,        iVECTOR   &dest  );         \
friend std::ostream &operator<<(std::ostream &output, const lVECTOR   &source);         \
friend std::istream &operator>>(std::istream &input,        lVECTOR   &dest  );

#endif
