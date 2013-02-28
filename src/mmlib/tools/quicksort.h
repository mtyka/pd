// PD is a free, modular C++ library for biomolecular simulation with a 
// flexible and scriptable Python interface. 
// Copyright (C) 2003-2013 Mike Tyka and Jon Rea
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __QUICKSORT_H
#define __QUICKSORT_H

#define LIMIT 5

template < class T >
void swap(T s[], int ind[], int i, int j){
	T tmp;
	int tmp2;

	tmp = s[i];
	s[i] = s[j];
	s[j] = tmp;

	tmp2 = ind[i];
	ind[i] = ind[j];
	ind[j] = tmp2;
}

template < class T >
void Insertionsort(T s[], int ind[], int length){
	int i, j;
	T elem;
	int index;

	for(i = 1; i < length; i++) {
		elem = s[i];
		index = ind[i];
		j = i;
		for(; j > 0 && elem < s[j - 1]; j--) {
			s[j] = s[j - 1];
			ind[j] = ind[j - 1];
		}
		s[j] = elem;
		ind[j] = index;
	}
}


template < class T >
void quicksort(T s[], int ind[], int left, int right){
	int i, j;
	T pivot;
	int middle;

	// switch to Insertionsort if the array to be sorted is
	// small enough.

	if(left + LIMIT > right)
		Insertionsort(s + left, ind + left, right - left + 1);
	else {
		middle = (left + right) / 2;
		if(s[middle] < s[left])
			swap(s, ind, left, middle);
		if(s[right] < s[left])
			swap(s, ind, right, left);
		if(s[right] < s[middle])
			swap(s, ind, right, middle);

		// The pivot will now be placed in position right-1

		swap(s, ind, right - 1, middle);

		pivot = s[right - 1];
		i = left;
		j = right - 1;

		for(;;) {
			while(s[++i] < pivot) {
			} // Don't modify ever!
			while(s[--j] > pivot) {
			} // Don't modify ever!

			if(i < j)
				swap(s, ind, i, j);
			else
				break;
		}

		// Restore the pivot

		swap(s, ind, right - 1, i);

		quicksort(s, ind, left, i - 1);
		quicksort(s, ind, i + 1, right);
	}
}

template < class T >
void qcksort(T s[], int ind[], int length){
	quicksort(s, ind, 0, length - 1);
}


template < class T >
inline void swap(T s[], int i, int j){
	T tmp;
	tmp = s[i];
	s[i] = s[j];
	s[j] = tmp;
}

template < class T >
inline void Insertionsort(T s[], int length){
	int i, j;
	T elem;

	for(i = 1; i < length; i++) {
		elem = s[i];
		j = i;
		for(; j > 0 && elem < s[j - 1]; j--) {
			s[j] = s[j - 1];
		}
		s[j] = elem;
	}
}


template < class T >
void quicksort(T s[], int left, int right){
	int i, j;
	T pivot;
	int middle;

	// switch to Insertionsort if the array to be sorted is
	// small enough.

	if(left + LIMIT > right)
		Insertionsort(s + left, right - left + 1);
	else {
		middle = (left + right) / 2;
		if(s[middle] < s[left])
			swap(s, left, middle);
		if(s[right] < s[left])
			swap(s, right, left);
		if(s[right] < s[middle])
			swap(s, right, middle);

		// The pivot will now be placed in position right-1
		swap(s, right - 1, middle);

		pivot = s[right - 1];
		i = left;
		j = right - 1;

		for(;;) {
			while(s[++i] < pivot) {} // Don't modify ever!
			while(s[--j] > pivot) {} // Don't modify ever!

			if(i < j)				swap(s, i, j);
			else    				break;
		}

		// Restore the pivot

		swap(s, right - 1, i);

		quicksort(s, left, i - 1);
		quicksort(s, i + 1, right);
	}
}

template < class T >
void qcksort(T s[], int length){
	quicksort(s, 0, length - 1);
}

#endif

