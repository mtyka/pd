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

#include "global.h"

#include "quicksort.h"

#define LIMIT 5

void swap(double s[], int ind[], int i, int j){
	double tmp;
	int tmp2;

	tmp = s[i];
	s[i] = s[j];
	s[j] = tmp;

	tmp2 = ind[i];
	ind[i] = ind[j];
	ind[j] = tmp2;
}

void Insertionsort(double s[], int ind[], int length){
	int i, j;
	double elem;
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


void quicksort(double s[], int ind[], int left, int right){
	int i, j;
	double pivot;
	int middle;

	// Switch to Insertionsort if the array to be sorted is
	// small enough.

	if(left + LIMIT > right)
		Insertionsort(s + left, ind + left, right - left + 1);
	else {

		// The pivot will be picked using the Median-of-Three strategy.

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

void qcksort(double s[], int ind[], int length){
	quicksort(s, ind, 0, length - 1);
}


// Integer versions of the above code


void iswap(int *s, int *ind, int i, int j){
	int tmp;
	int tmp2;

	tmp = s[i];
	s[i] = s[j];
	s[j] = tmp;

	tmp2 = ind[i];
	ind[i] = ind[j];
	ind[j] = tmp2;
}

void iInsertionsort(int *s, int *ind, int length){
	int i, j;
	int elem;
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


void iquicksort(int *s, int *ind, int left, int right){
	int i, j;
	int pivot;
	int middle;

	// Switch to Insertionsort if the array to be sorted is
	// small enough.

	if(left + LIMIT > right)
		iInsertionsort(s + left, ind + left, right - left + 1);
	else {

		// The pivot will be picked using the Median-of-Three strategy.

		middle = (left + right) / 2;
		if(s[middle] < s[left])
			iswap(s, ind, left, middle);
		if(s[right] < s[left])
			iswap(s, ind, right, left);
		if(s[right] < s[middle])
			iswap(s, ind, right, middle);

		// The pivot will now be placed in position right-1

		iswap(s, ind, right - 1, middle);

		pivot = s[right - 1];
		i = left;
		j = right - 1;

		for(;;) {
			while(s[++i] < pivot) {
			} // Don't modify ever!
			while(s[--j] > pivot) {
			} // Don't modify ever!

			if(i < j)
				iswap(s, ind, i, j);
			else
				break;
		}

		// Restore the pivot

		iswap(s, ind, right - 1, i);

		iquicksort(s, ind, left, i - 1);
		iquicksort(s, ind, i + 1, right);
	}
}

void iqcksort(int *s, int *ind, int length){
	iquicksort(s, ind, 0, length - 1);
}
