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

#ifndef __RMS_H
#define __RMS_H

// required in headers as template bodies must be in the header
#include "maths/maths.h"
#include "tools/quicksort.h" 

/// This is a freestanding function calculating the cRMS between two supplied Maths::Tvector<T> arrays
/// It calculates the cRMS between two structures without need for iterative superimposition
/// according to an algorithm developed by Kabsch and uses two supplied arrays
/// as the coordinates (vsource & vtarget)
///
/// W. KABSCH, Acta Cryst. (1978). A34, 827-828
/// A discussion o f the solution for the best rotation to relate two sets of vectors.
/// See also KABSCH, W. (1976). Acta Cryst. A32, 922-923. for proofs.
///
template < class T >
T calcVectorCRMS_superimp(Maths::Tvector < T > *vsource, const Maths::Tvector < T > *vtarget, int nvectors){
	Maths::Tvector < T > source_cog;
	Maths::Tvector < T > tcog;
	Maths::matrix3x3 rmat;
	T cRMS = calcVectorCRMS(vsource, vtarget, nvectors, &source_cog, &tcog, &rmat);
	for(int i = 0; i < nvectors; i++) { // change the source structure
		vsource[i].sub(source_cog);
		vsource[i].mulmat(rmat);
		vsource[i].add(tcog);
	}
	return cRMS;
}

template < class T >
T calcVectorCRMS(const Maths::Tvector < T > *vsource, const Maths::Tvector < T > *vtarget, int nvectors, Maths::Tvector < T > *source_cog, Maths::Tvector < T > *target_cog, Maths::matrix3x3 * endU){
	if(nvectors <= 0)
		return -1.0;
	int i;
	Maths::Tvector < T > Tvec; // translation to be applied to atom[i].p to remove translational differences
	Maths::Tvector < T > p(0, 0, 0), t(0, 0, 0);
	Maths::Tvector < T > pcog(0, 0, 0), tcog(0, 0, 0); // centre of geometries

	Maths::matrix3x3 R, Rt;
	T E0;
	T crms = 0;

	Maths::dvector v[3];
	double l1, l2, l3;
	T lsort[3];
	int lindex[3] = { 0, 1, 2 };
	Maths::dvector a1, a2, a3;
	Maths::dvector b1, b2, b3;
	Maths::matrix3x3 U;

	// find the translation Maths::Tvector<T>
	for(i = 0; i < nvectors; i++) {
		pcog.add(vsource[i]);
		tcog.add(vtarget[i]);
	}
	tcog.div((T) nvectors);
	pcog.div((T) nvectors);
	Tvec.diff(tcog, pcog);

	// find initial E0 and R
	E0 = 0;
	R.setToNull();
	for(i = 0; i < nvectors; i++) {
		p.setTo(vsource[i]);
		t.setTo(vtarget[i]);

		E0 += sqr((p.x - pcog.x) - (t.x - tcog.x)) +
			sqr((p.y - pcog.y) - (t.y - tcog.y)) + sqr((p.z - pcog.z) - (t.z - tcog.z));

		R.r[0][0] += (t.x - tcog.x) * (p.x - pcog.x);
		R.r[1][0] += (t.y - tcog.y) * (p.x - pcog.x);
		R.r[2][0] += (t.z - tcog.z) * (p.x - pcog.x);
		R.r[0][1] += (t.x - tcog.x) * (p.y - pcog.y);
		R.r[1][1] += (t.y - tcog.y) * (p.y - pcog.y);
		R.r[2][1] += (t.z - tcog.z) * (p.y - pcog.y);
		R.r[0][2] += (t.x - tcog.x) * (p.z - pcog.z);
		R.r[1][2] += (t.y - tcog.y) * (p.z - pcog.z);
		R.r[2][2] += (t.z - tcog.z) * (p.z - pcog.z);
	}

	if(E0 < 0.001)
		return sqrt(E0 / (T) nvectors); // dont do a rotation if the cRMS is already so low

	Rt.setTo(R);
	Rt.transpose();
	Rt.postmul(R);
	Rt.diagonaliseSymetric(l1, l2, l3, v[0], v[1], v[2]);

	// sort
	lsort[0] = (T) l1;
	lsort[1] = (T) l2;
	lsort[2] = (T) l3;
	qcksort(lsort, lindex, 3);

	a1.setTo(v[lindex[2]]); // set a1,a2 according to l1>l2>l3
	a2.setTo(v[lindex[1]]);
	a3.crossProduct(a1, a2);
	a1.unify();
	a2.unify();
	a3.unify();

	b1.setTo(a1);
	b1.mulmat(R);
	b1.unify();
	b2.setTo(a2);
	b2.mulmat(R);
	b2.unify();
	b3.crossProduct(b1, b2);

	U.r[0][0] = b1.x * a1.x + b2.x * a2.x + b3.x * a3.x;
	U.r[0][1] = b1.x * a1.y + b2.x * a2.y + b3.x * a3.y;
	U.r[0][2] = b1.x * a1.z + b2.x * a2.z + b3.x * a3.z;

	U.r[1][0] = b1.y * a1.x + b2.y * a2.x + b3.y * a3.x;
	U.r[1][1] = b1.y * a1.y + b2.y * a2.y + b3.y * a3.y;
	U.r[1][2] = b1.y * a1.z + b2.y * a2.z + b3.y * a3.z;

	U.r[2][0] = b1.z * a1.x + b2.z * a2.x + b3.z * a3.x;
	U.r[2][1] = b1.z * a1.y + b2.z * a2.y + b3.z * a3.y;
	U.r[2][2] = b1.z * a1.z + b2.z * a2.z + b3.z * a3.z;

	// find numerical cRMS
	crms = 0.0;
	for(i = 0; i < nvectors; i++) {
		p.setTo(vsource[i]); // load atom
		p.sub(pcog); // translate COG to origin
		p.mulmat(U);
		t.setTo(vtarget[i]);
		crms += sqr(p.x - t.x + tcog.x) + sqr(p.y - t.y + tcog.y) + sqr(p.z - t.z + tcog.z);
	}

	crms = sqrt(crms / (T) nvectors);

	if(source_cog != NULL)
		source_cog->setTo(pcog);
	if(target_cog != NULL)
		target_cog->setTo(tcog);
	if(endU != NULL)
		endU->setTo(U);

	return crms;
}


template < class T >
T calcVectorCRMS_ignoreInvalid(const Maths::Tvector < T > *vsource, const Maths::Tvector < T > *vtarget, int nvectors, Maths::Tvector < T > *source_cog, Maths::Tvector < T > *target_cog, Maths::matrix3x3 * endU)
{
	if(nvectors <= 0)
		return -1.0;
	int i;
	Maths::Tvector < T > Tvec; // translation to be applied to atom[i].p to remove translational differences
	Maths::Tvector < T > p(0, 0, 0), t(0, 0, 0);
	Maths::Tvector < T > pcog(0, 0, 0), tcog(0, 0, 0); // centre of geometries

	Maths::matrix3x3 R, Rt;
	T E0;
	T crms = 0;

	Maths::dvector v[3];
	double l1, l2, l3;
	T lsort[3];
	int lindex[3] = { 0, 1, 2 };
	Maths::dvector a1, a2, a3;
	Maths::dvector b1, b2, b3;
	Maths::matrix3x3 U;

	int nvalidvectors = 0;

	// find the translation Maths::Tvector<T>
	for(i = 0; i < nvectors; i++) {
		if(vsource[i].magnitude < 0)
			continue;
		if(vtarget[i].magnitude < 0)
			continue;
		pcog.add(vsource[i]);
		tcog.add(vtarget[i]);
		nvalidvectors++;
	}
	if(nvalidvectors <= 0)
		return -1.0;
	tcog.div((T) nvalidvectors);
	pcog.div((T) nvalidvectors);
	Tvec.diff(&tcog, &pcog);

	// find initial E0 and R
	E0 = 0;
	R.setToNull();
	for(i = 0; i < nvectors; i++) {
		if(vsource[i].magnitude < 0)
			continue;
		if(vtarget[i].magnitude < 0)
			continue;
		p.setTo(vsource[i]);
		t.setTo(vtarget[i]);

		E0 += sqr((p.x - pcog.x) - (t.x - tcog.x)) +
			sqr((p.y - pcog.y) - (t.y - tcog.y)) + sqr((p.z - pcog.z) - (t.z - tcog.z));

		R.r[0][0] += (t.x - tcog.x) * (p.x - pcog.x);
		R.r[1][0] += (t.y - tcog.y) * (p.x - pcog.x);
		R.r[2][0] += (t.z - tcog.z) * (p.x - pcog.x);
		R.r[0][1] += (t.x - tcog.x) * (p.y - pcog.y);
		R.r[1][1] += (t.y - tcog.y) * (p.y - pcog.y);
		R.r[2][1] += (t.z - tcog.z) * (p.y - pcog.y);
		R.r[0][2] += (t.x - tcog.x) * (p.z - pcog.z);
		R.r[1][2] += (t.y - tcog.y) * (p.z - pcog.z);
		R.r[2][2] += (t.z - tcog.z) * (p.z - pcog.z);
	}

	if(E0 < 0.001)
		return sqrt(E0 / (T) nvalidvectors); // dont do a rotation if the cRMS is already so low

	Rt.setTo(R);
	Rt.transpose();
	Rt.postmul(R);
	Rt.diagonaliseSymetric(l1, l2, l3, v[0], v[1], v[2]);

	// sort
	lsort[0] = (T) l1;
	lsort[1] = (T) l2;
	lsort[2] = (T) l3;
	qcksort(lsort, lindex, 3);

	a1.setTo(v[lindex[2]]); // set a1,a2 according to l1>l2>l3
	a2.setTo(v[lindex[1]]);
	a3.crossProduct(a1, a2);
	a1.unify();
	a2.unify();
	a3.unify();

	b1.setTo(a1);
	b1.mulmat(R);
	b1.unify();
	b2.setTo(a2);
	b2.mulmat(R);
	b2.unify();
	b3.crossProduct(b1, b2);

	U.r[0][0] = b1.x * a1.x + b2.x * a2.x + b3.x * a3.x;
	U.r[0][1] = b1.x * a1.y + b2.x * a2.y + b3.x * a3.y;
	U.r[0][2] = b1.x * a1.z + b2.x * a2.z + b3.x * a3.z;

	U.r[1][0] = b1.y * a1.x + b2.y * a2.x + b3.y * a3.x;
	U.r[1][1] = b1.y * a1.y + b2.y * a2.y + b3.y * a3.y;
	U.r[1][2] = b1.y * a1.z + b2.y * a2.z + b3.y * a3.z;

	U.r[2][0] = b1.z * a1.x + b2.z * a2.x + b3.z * a3.x;
	U.r[2][1] = b1.z * a1.y + b2.z * a2.y + b3.z * a3.y;
	U.r[2][2] = b1.z * a1.z + b2.z * a2.z + b3.z * a3.z;

	// find numerical cRMS
	crms = 0.0;
	for(i = 0; i < nvectors; i++) {
		if(vsource[i].magnitude < 0)
			continue;
		if(vtarget[i].magnitude < 0)
			continue;
		p.setTo(vsource[i]); // load atom
		p.sub(pcog); // translate COG to origin
		p.mulmat(U);
		t.setTo(vtarget[i]);
		crms += sqr(p.x - t.x + tcog.x) + sqr(p.y - t.y + tcog.y) + sqr(p.z - t.z + tcog.z);
	}

	crms = sqrt(crms / (T) nvalidvectors);

	if(source_cog != NULL)
		source_cog->setTo(pcog);
	if(target_cog != NULL)
		target_cog->setTo(tcog);
	if(endU != NULL)
		endU->setTo(U);

	return crms;
}


template < class T >
T calcVectorCRMS(const Maths::Tvector < T > *vsource, const Maths::Tvector < T > *vtarget, int nvectors){
	return calcVectorCRMS(vsource, vtarget, nvectors, (Maths::Tvector < T > *)NULL, (Maths::Tvector < T > *)NULL, NULL);
}


// This is a freestanding function calculating the dRMS between two supplied Maths::Tvector<T> arrays
template < class T >
T calcVectorDRMS(Maths::Tvector < T > *vsource, Maths::Tvector < T > *vtarget, int nvectors){
	T drms = 0.0;
	int i, j;
	int count = 0;
	T dista, distb;

	for(i = 0; i < nvectors; i++) { // calculate dRMS
		for(j = i + 1; j < nvectors; j++) {
			dista = vsource[i].dist(vsource[j]);
			distb = vtarget[i].dist(vtarget[j]);
			drms += sqr(dista - distb);
			count++;
		}
	}
	drms /= (T) count;
	drms = sqrt(drms);

	return drms;
}


#endif

