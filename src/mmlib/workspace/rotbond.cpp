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

#include "rotbond.h"

#include "forcefields/ffparam.h"
#include "workspace/workspace.h"

using namespace Maths;

// Small helper function for std::qsort used in 'ReAllocate'
int intcmp(const void *elem1, const void *elem2){
	if((*(int *) elem1) < (*(int *) elem2))
		return -1;
	if((*(int *) elem1) > (*(int *) elem2))
		return 1;
	return 0;
}

RotBond::RotBond():
	WorkSpaceComponentBase(NULL)
{
	m_OutputLevel = Verbosity::Normal;
}

void RotBond::clear()
{
	rotbond.clear();
}

size_t RotBond::memuse(int level)
{
	return sizeof(*this) + rotbond.size()*sizeof(RotatableBond);
}

void RotBond::reinit( WorkSpace* _wspace )
{
	if(_wspace==NULL) THROW(ArgumentException,"RotBond::ReAllocate() given _wspace==NULL");
	wspace = _wspace;

	if( m_OutputLevel ) printf("Assembling rotatable bond list ... ");

	rotbond.clear(); // clear the list

	// Some proxy variables
	ParticleStore  &atom = wspace->atom;
	const size_t natom = wspace->atom.size();
	const size_t nmolecules = wspace->mol.size();

	int bi;
	int bj;
	unsigned i, j;
	int iat, n, c;
	int *list = new int[natom + 1];
	int* localStatus = new int[natom];
	int listentries = 1;
	int curentry = -1;
	int listatom;
	int neighatom;
	int natomstooneside;
	int sideofbond;
	bool bondinring;

	for(iat = 0; iat < natom; iat++) {
		for(c = 0; c < atom[iat].cov12atom.size(); c++) { //collect all the covalent links, assume they're all bond...

			if(iat < atom[iat].cov12atom[c].i)
				continue; //only record half (i.e. don't double record 1-3 and 3-1)

			RotatableBond newrotbond;
			newrotbond.i = iat; //record atoms
			newrotbond.j = atom[iat].cov12atom[c].i;

			for(sideofbond = 0; sideofbond < 2; sideofbond++) { //try either way (to chose the smaller side)
				natomstooneside = 0;
				bondinring = false;

				if(sideofbond == 0) {
					bi = newrotbond.i;
					bj = newrotbond.j;
				} else {
					bj = newrotbond.i;
					bi = newrotbond.j;
				}

				// count how many atoms 'hang' on either side of the bond
				// and see if it's a bond in a ring (in which case it's non rotatable)

				for(i = 0; i < natom; i++)
					localStatus[i] = -1;// mark as unrotated
				
				localStatus[bj] = 1; // mark i as rotated (so we cant traverse the i-j bond

				list[0] = bj; // set j as the first entry
				listentries = 1;
				curentry = -1;

				while(curentry < (listentries - 1)) {
					curentry++;
					listatom = list[curentry];
					if(localStatus[listatom] > 1)
						continue; // ignore if already rotated
					localStatus[listatom] = 2;// mark as dealt with
					// now find listatom's covalently linked atoms to create new entries to the list
					int atom_listatom_cov12atom_size = atom[listatom].cov12atom.size();
					for(n = 0; n < atom_listatom_cov12atom_size; n++) {
						neighatom = atom[listatom].cov12atom[n].i; // get atom index of cov neighbour
						if(neighatom == listatom)
							continue; // ignore direct backstepping
						if(localStatus[neighatom] > 0)
							continue; // ignore if already rotated
						if(listatom == bj) {
							if(neighatom == bi)
								continue; // do not Step through the bond onto the other side
						} else {
							if(neighatom == bi) {
								bondinring = true;
								break;
							} // do not Step through the bond onto the other side
						}

						localStatus[neighatom] = 1;// mark as listed
						// if arrive here - add to list for processing
						list[listentries] = neighatom;
						// printf("added atom %d \n",neighatom);
						listentries++;
						if(listentries >= (natom + 1)) {
							printf(" Error!!\n");
							THROW(ProcedureException,"List overflow in RotBond::rotateBond(...) \n");
						}
					}
					if(bondinring)
						break;
				}
				// if the number of listentries is smaller than half the atoms in this
				// molecule then finish - we've found the smaller side of the bond
				if(listentries < (natom / 2))
					break;
			}

			// now analyse list and find the segments of atoms to be
			// the atoms to be rotated

			// ignore bonds that are part of rings
			if(bondinring) continue;

			//no point in rotating if there's only one atom there anyway
			if(listentries <= 1) continue;

			qsort((void *) &list[0], listentries, sizeof(int), intcmp);

			int segments = 0;

			newrotbond.segmentstart[segments] = list[0];
			for(i = 1; i < listentries; i++) { // going through all the selected atoms
				if(list[i] == (list[i - 1] + 1))
					continue;// no atom number gap so continue;
				//ok end of segment
				newrotbond.segmentend[segments] = list[i - 1]; // finish off current segment
				segments++;
				if(segments >= (MAXROTATABLEBONDSEGS - 1)) {
					delete[]list; // 1st free list memory
					THROW(CodeException,"MAXROTATABLEBONDSEGS Overflow in ReAllocate()\n");
				}
				newrotbond.segmentstart[segments] = list[i]; // and start a new one
			}
			newrotbond.segmentend[segments] = list[i - 1]; // finish off current segment
			segments++;

			const char *itype = atom[bi].pdbname.c_str();;
			const char *jtype = atom[bj].pdbname.c_str();;

			// determine Type of bond

			if(((strcmp(itype, "N") == 0) && (strcmp(jtype, "C") == 0)) || ((strcmp(itype, "C") == 0) && (strcmp(jtype, "N") == 0))) { //printf("Peptide Bond, Planar\n");
				newrotbond.Type = RotatableBond::Omega;
			}

			if(((strcmp(itype, "N") == 0) && (strcmp(jtype, "CA") == 0)) || ((strcmp(itype, "CA") == 0) && (strcmp(jtype, "N") == 0))) { //printf("Phi\n");
				newrotbond.Type = RotatableBond::Phi;
			}

			if(((strcmp(itype, "C") == 0) && (strcmp(jtype, "CA") == 0)) || ((strcmp(itype, "CA") == 0) && (strcmp(jtype, "C") == 0))) { //printf("Psi\n");
				newrotbond.Type = RotatableBond::Psi;
			}

			// chose two arbitary pi,pj atoms:
			int ni, nj;
			newrotbond.ip = -1;
			for(ni = 0; ni < atom[newrotbond.i].cov12atom.size(); ni++) {
				if(atom[newrotbond.i].cov12atom[ni].i != newrotbond.j) {
					newrotbond.ip = atom[newrotbond.i].cov12atom[ni].i;
					break;
				}
			}
			newrotbond.jp = -1;
			for(nj = 0; nj < atom[newrotbond.j].cov12atom.size(); nj++) {
				if(atom[newrotbond.j].cov12atom[nj].i != newrotbond.i) {
					newrotbond.jp = atom[newrotbond.j].cov12atom[nj].i;
					break;
				}
			}

			rotbond.push_back(newrotbond);
			//printf("%3d %3d %3d %3d \n",newrotbond.ip,newrotbond.i,newrotbond.j,newrotbond.jp );

			calcTorsion(rotbond.size()-1);
		}
	}

	delete[]localStatus;
	delete[]list; // free list memory

	// 'connect' phi & psi for concerted changes
	for(i = 0; i < rotbond.size(); i++) {
		if(rotbond[i].Type == RotatableBond::Psi) {
			for(j = 0; j < rotbond.size(); j++) {
				if((rotbond[j].Type == RotatableBond::Phi) && // find phi
					(atom[rotbond[i].i].ir == atom[rotbond[j].i].ir)) // of the same residue
					break;
			}
			if(j < rotbond.size())
				rotbond[i].status = j; // remember the associated phi angle
			else
				rotbond[i].status = -1; // cant find associated phi
		};
		if(rotbond[i].Type == RotatableBond::Phi) {
			for(j = 0; j < rotbond.size(); j++) {
				if((rotbond[j].Type == RotatableBond::Psi) && // find phi
					(atom[rotbond[i].i].ir == atom[rotbond[j].i].ir)) // of the same residue
					break;
			}
			if(j < rotbond.size())
				rotbond[i].status = j; // remember the associated psi angle
			else
				rotbond[i].status = -1; // cant find associated psi
		}
	}

	if( m_OutputLevel ) 
	{	
		printf(" Done!\n");
		printf("Internal degrees of Freedom:\n");
		printf("  Cartesian: %d\n", Maths::max((int)0,(int)(natom * 3 - 6)));
		printf("  Torsional: %d\n", Maths::max((int)0,(int)(rotbond.size() + (nmolecules - 1) * 6)));
	}
}

int RotBond::find(int atomI, int atomJ )
{
	for( size_t q = 0; q < rotbond.size(); q++ )
	{
		if( ( rotbond[q].i == atomI ) && ( rotbond[q].j == atomJ ) )
		{
			return q;
		}
		else if( ( rotbond[q].i == atomJ ) && ( rotbond[q].j == atomI ) )
		{
			return q;
		}
	}
	return -1;
}

int RotBond::rotate(int irbond, double angle){

	int natom = wspace->atom.size();
	SnapShotAtom* atom = wspace->cur.atom;

	int i, s;
	matrix3x3 rmat;
	dvector origin(atom[rotbond[irbond].i].p);
	dvector axis;

	axis.setTo(atom[rotbond[irbond].j].p);
	axis.sub(atom[rotbond[irbond].i].p);
	rmat.setToAxisRot(axis, angle);

	s = 0;
	while(rotbond[irbond].segmentstart[s] >= 0) {
		for(i = rotbond[irbond].segmentstart[s]; i <= rotbond[irbond].segmentend[s]; i++) {
			atom[i].p.sub(origin); // translate to axis origin
			atom[i].p.mulmat(rmat); // rotate
			atom[i].p.add(origin); // translate back to original position
		}
		s++;
		if(s >= MAXROTATABLEBONDSEGS)
			break;
	}

	return 0;
}

int RotBond::rotate(int irbond, double angle, int maxresidue){

	int natom = wspace->atom.size();
	SnapShotAtom* atom = wspace->cur.atom;

	int i, s;
	int bondresidue;

	matrix3x3 rmat;
	dvector origin(atom[rotbond[irbond].i].p);
	dvector axis;

	bondresidue = wspace->atom[rotbond[irbond].j].ir;

	axis.setTo(atom[rotbond[irbond].j].p);
	axis.sub(atom[rotbond[irbond].i].p);
	rmat.setToAxisRot(axis, angle);

	s = 0;
	while(rotbond[irbond].segmentstart[s] >= 0) {
		for(i = rotbond[irbond].segmentstart[s]; i <= rotbond[irbond].segmentend[s]; i++) {
			if(abs(wspace->atom[i].ir - bondresidue) > maxresidue)
				continue;
			atom[i].p.sub(origin); // translate to axis origin
			atom[i].p.mulmat(rmat); // rotate
			atom[i].p.add(origin); // translate back to original position
		}
		s++;
		if(s >= MAXROTATABLEBONDSEGS)
			break;
	}

	return 0;
}

double RotBond::calcTorsion(int irbond){

	if(rotbond[irbond].ip < 0)
		return 0.0;
	if(rotbond[irbond].i < 0)
		return 0.0;
	if(rotbond[irbond].j < 0)
		return 0.0;
	if(rotbond[irbond].jp < 0)
		return 0.0;


	// Some proxy variables
	const ParticleStore  &atom = wspace->atom;
	size_t natom = wspace->atom.size();
	size_t nmolecules = wspace->mol.size();

	rotbond[irbond].lastphi = rotbond[irbond].phi;
	rotbond[irbond].phi = (double) calcTorsionAngle(
		wspace->cur.atom[rotbond[irbond].ip].p, 
		wspace->cur.atom[rotbond[irbond].i].p,
		wspace->cur.atom[rotbond[irbond].j].p, 
		wspace->cur.atom[rotbond[irbond].jp].p);

	return rotbond[irbond].phi;
}

void RotBond::silence( Verbosity::Type _OutputLevel ) 
{
	m_OutputLevel = _OutputLevel; 
}

bool RotBond::isDelocalised( int rotBondIndex )
{
	// Normally 'rotatable' single bonds also include those that are delocalised,
	// i.e. the C-N bonds at the end of a glutamine or agrinine. These should not nececarily be included in a torsional model.
	// The way that I have decided to isolate these bonds is by looking at the 'IMPROPER' definitions in the forcefield file.
	// This should define those bonds that are currently marked as 'single' but should not nececarily be interpreted as such.

	int natom = wspace->atom.size();
	int nresidues = wspace->res.size();
	const ParticleStore& atom = wspace->atom;
	const ResidueStore& res = wspace->res;

	const MoleculeDefinition *molDef;
	RotatableBond &rotBondRef = rotbond[rotBondIndex];

	for( int i = 0; i < nresidues; i++ )
	{
		int psOffset = res[i].ifirst;

		if( rotBondRef.i >= res[i].ifirst && rotBondRef.i <= res[i].ilast &&
			rotBondRef.j >= res[i].ifirst && rotBondRef.j <= res[i].ilast )
		{
			// the rotbond is in this pResidue
			molDef = res[i].param;
			int hasBothCount = 0; // the number of improper definitions containing both of the atoms.
			for( unsigned j = 0; j < molDef->improper.size(); j++ )
			{
				if( rotBondRef.i == (molDef->findAtomRaw( molDef->improper[j].ani ) + psOffset) )
				{
					if ( rotBondRef.j == (molDef->findAtomRaw( molDef->improper[j].anj ) + psOffset) ) hasBothCount++;
					else if( rotBondRef.j == (molDef->findAtomRaw( molDef->improper[j].ana ) + psOffset) ) hasBothCount++;
					else if( rotBondRef.j == (molDef->findAtomRaw( molDef->improper[j].anb ) + psOffset) ) hasBothCount++;
				}
				else if( rotBondRef.i == (molDef->findAtomRaw( molDef->improper[j].anj ) + psOffset) )
				{
					if ( rotBondRef.j == (molDef->findAtomRaw( molDef->improper[j].ani ) + psOffset) ) hasBothCount++;
					else if( rotBondRef.j == (molDef->findAtomRaw( molDef->improper[j].ana ) + psOffset) ) hasBothCount++;
					else if( rotBondRef.j == (molDef->findAtomRaw( molDef->improper[j].anb ) + psOffset) ) hasBothCount++;
				}
				else if( rotBondRef.i == (molDef->findAtomRaw( molDef->improper[j].ana ) + psOffset) )
				{
					if ( rotBondRef.j == (molDef->findAtomRaw( molDef->improper[j].ani ) + psOffset) ) hasBothCount++;
					else if( rotBondRef.j == (molDef->findAtomRaw( molDef->improper[j].anj ) + psOffset) ) hasBothCount++;
					else if( rotBondRef.j == (molDef->findAtomRaw( molDef->improper[j].anb ) + psOffset) ) hasBothCount++;
				}
				else if( rotBondRef.i == (molDef->findAtomRaw( molDef->improper[j].anb ) + psOffset) )
				{
					if ( rotBondRef.j == (molDef->findAtomRaw( molDef->improper[j].ani ) + psOffset) ) hasBothCount++;
					else if( rotBondRef.j == (molDef->findAtomRaw( molDef->improper[j].ana ) + psOffset) ) hasBothCount++;
					else if( rotBondRef.j == (molDef->findAtomRaw( molDef->improper[j].anj ) + psOffset) ) hasBothCount++;
				}

				if( hasBothCount == 2 )
					return true; // we have found two matching improper definitions. We therefore assume that this bond is delocalised.
			}
			return false; // less than two impropers including these atoms were found.
		}
	}
	THROW(CodeException,"Assumption error in isDelocalised(). Atom index seek was unsucessful!");
}

// -----------------------------------------------------------------------------------------------------------
// Begin: 'RotationDefinition'
// -----------------------------------------------------------------------------------------------------------

void RotationDefinition::performRotation( double desiredAngle )
{
	D_ASSERT( isValid(), CodeException, "RotationDefinition is not valid!");

	matrix3x3 rmat;
	dvector axis;
	dvector origin( *atom2 );

	double rotateCurrentAngle = getCurrentTorsionAngle();
	rotateCurrentAngle = desiredAngle - rotateCurrentAngle; // and now ... what do we need to rotate by?

	// Check that we are actually making a change!
	if( Maths::SigFigEquality( rotateCurrentAngle, 0.0, 8 ) )
	{
		return;
	}

	// assign the rotation matrix
	axis.setTo(*atom3);
	axis.sub(*atom2);
	rmat.setToAxisRot( axis, rotateCurrentAngle );

	// rotate all the atoms that are affetced by the rotation within the Scope of the loop definition
	for( int i = startAtomIndex; i <= endAtomIndex; i++)
	{
		molBase->atomxyz(i).sub(origin); // translate to axis origin
		molBase->atomxyz(i).mulmat(rmat); // rotate position
		molBase->atomxyz(i).add(origin); // translate back to original position
	}

	if( arseAtom > 0 )
	{
		molBase->atomxyz(arseAtom).sub(origin); // translate to axis origin
		molBase->atomxyz(arseAtom).mulmat(rmat); // rotate position
		molBase->atomxyz(arseAtom).add(origin); // translate back to original position
	}
}

void RotationDefinition::perturbRotation( double deltaAngle )
{
	matrix3x3 rmat;
	dvector axis;
	dvector origin( *atom2 );

	// assign the rotation matrix
	axis.setTo(*atom3);
	axis.sub (*atom2);
	rmat.setToAxisRot( axis, deltaAngle );

	// rotate all the atoms that are affetced by the rotation within the Scope of the loop definition
	for( int i = startAtomIndex; i <= endAtomIndex; i++)
	{
		molBase->atomxyz(i).sub(origin); // translate to axis origin
		molBase->atomxyz(i).mulmat(rmat); // rotate position
		molBase->atomxyz(i).add(origin); // translate back to original position
	}

	if( arseAtom != -1 )
	{
		molBase->atomxyz(arseAtom).sub(origin); // translate to axis origin
		molBase->atomxyz(arseAtom).mulmat(rmat); // rotate position
		molBase->atomxyz(arseAtom).add(origin); // translate back to original position
	}

}

// -----------------------------------------------------------------------------------------------------------
// End: 'RotationDefinition'
// -----------------------------------------------------------------------------------------------------------



