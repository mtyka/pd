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

// Global pre-compiled header
#include "global.h"

// Self Header
#include "system/molecule.h"

// mmLib Includes
#include "fileio/pdb.h"
#include "forcefields/forcefield.h"

// Namespace Includes
using namespace Maths;

////////////////////////////////////////////////////////////////////////////////
//
// Inspectors and info functions


void MoleculeBase::info() const
{
	printf(" -------------------------------------\n");
	printf(" Atoms:        %8d\n", nAtoms());
	printf(" Residues:     %8d\n", nResidues());
	printf(" Molecules:    %8d\n", nMolecules());
	printf(" Total mass:   %8.4lf amu\n", getTotalMass() / Physics::PhysicsConst::amu);
	printf(" Total charge: %+8.2lf e\n", getTotalCharge() );
}

void MoleculeBase::detail( ) const 
{
	detail( PickAllParticles() );
}

void MoleculeBase::detail( const PickBase &picker ) const
{

	printf(".-------- Atom --------..-- Name --..- Parent -..----- Properties ------..-- Covalent --- \n");
 	printf("    no   res   mol  grp   pdb    ff  pdb   name   Z  mass  rad  chg  flags(HBCSK) \n");

	for(size_t i = 0; i < nAtoms(); i++) 
	{
		if( !picker.matches( atom[i] ) ) continue;

		printf("%6d%6d%6d%6d  %4s%6s %3s %6s%3d %6.3lf %5.3lf % 7.4lf   ",
			atom[i].i, 
			atom[i].ir, 
			atom[i].imol,
			atom[i].igroup,
			atom[i].pdbname.c_str(), 
			atom[i].rawname.c_str(),
			atom[i].parentl3name.c_str(),
			atom[i].parentname.c_str(),
			atom[i].Z,
			atom[i].mass,
			atom[i].radius,
			atom[i].charge
		);

		// flags
		printf("%c", atom[i].isHydrogen()      ? 'H' : '-' );
		printf("%c", atom[i].isBackbone()      ? 'B' : '-' );
		printf("%c", atom[i].isCAlpha()        ? 'C' : '-' );
		printf("%c", atom[i].isStatic()        ? 'S' : '-' );
		printf("%c", atom[i].isKnownStructure()? 'K' : '-' );
		printf("  ");

		for(size_t r = 0; r < atom[i].cov12atom.size(); r++) 
		{
			printf(" %d(%s) ", atom[i].cov12atom[r].i, atom[atom[i].cov12atom[r].i].pdbname.c_str());
		}

		printf(" \n");
	}
}

void MoleculeBase::printBondConnectivity() const
{
	for(size_t iat=0;iat<nAtoms();iat++){
		printf("%d: ",iat);
		for(size_t k=0;k<atom[iat].cov12atom.size();k++){
			printf(" %d ",atom[iat].cov12atom[k].i);
		}
		printf("\n");
	}
}

int MoleculeBase::findParticleBy_ffname(int ir,const std::string &name) const 
{
	for(size_t i = 0; i < nAtoms(); i++) 
    {
		if(ir == atom[i].ir)
			if(strcmp(name.c_str(), atom[i].rawname.c_str()) == 0)
				return i;
	}
	return -1;
}

int MoleculeBase::findParticle(int ir,const std::string &name) const 
{
	for(size_t i = 0; i < nAtoms(); i++) {
		if(ir == atom[i].ir)
			if(strcmp(name.c_str(), atom[i].pdbname.c_str()) == 0)
				return i;
	}
	return -1;
}













////////////////////////////////////////////////////////////////////////////////
//
// Setup Functions 

int MoleculeBase::createCovalentStructure(){
	unsigned ir,iat;
	unsigned ibnd;
	int tindex,pindex;
	int errorstatus=0;
	for(ir = 0; ir < nResidues(); ir++) {
		// add bonds (i.e. covalent structure)

		// for every atom in the residue go through all the bonds and add them
		for(iat = 0; iat < res[ir].param->atom.size(); iat++) {
			// find the index of the parent particle
			tindex = findParticleBy_ffname(ir, res[ir].param->atom[iat].rawname.c_str());
			for(ibnd = 0; ibnd < res[ir].param->atom[iat].r_covalent.size(); ibnd++) {
				pindex = findParticleBy_ffname(ir + res[ir].param->atom[iat].r_covalent[ibnd].roi,
					res[ir].param->atom[iat].r_covalent[ibnd].ani.c_str());
				if(pindex < 0) {
					printf("\nERROR: Cannot find particle %d:%s - unable to make bond\n",
						ir + res[ir].param->atom[iat].r_covalent[ibnd].roi,
						res[ir].param->atom[iat].r_covalent[ibnd].ani.c_str());
					errorstatus = 1;
					continue;
				}
				// add covalent link to mother atom
				CovalentAtom newcov;
				newcov.i = pindex;

				atom[tindex].cov12atom.push_back(newcov);
			}
		}
	}
	if(errorstatus != 0) throw( ProcedureException("Errors or inconsistencies during the creation of 1:2 covalent structure. \nCheck the forcefield input files and/or input sequences. Ensure that capping \nresidue types are used at the termini. "));
	if( 0 != calc1314BondOrders() ) throw( ProcedureException("Errors or inconsistencies during the creation of 1:3 and 1:4 covalent structure. \nCheck the forcefield input files and/or input sequences. Ensure that capping \nresidue types are used at the termini. "));
	return 0;
}


int MoleculeBase::checkParameters(){
	// check loadParameters() has been called (at least once)
	if(!loadedparams){
		printf("Loaded params has not been called\n");
		return -1;
	}
	return 0;
}	

int MoleculeBase::checkCovalentStructure(){
	unsigned j, k;
	int resti, ok;
	unsigned l;
	int return_value = 0;

	for(j = 0; j < nAtoms() ; j++) { // for each atom
		for(k = 0; k < atom[j].cov12atom.size(); k++) { // for every of it's bonds

			// get the atom reference number for this restraint
			resti = atom[j].cov12atom[k].i;

			// Now check if the atom linked to, also links back to this atom.
			ok = 0; // (false)
			for(l = 0; l < atom[resti].cov12atom.size(); l++) {
				// is atom j part of atom resti's
				// covalent list ?
				if(atom[resti].cov12atom[l].i == j) { // Yes ?
					ok = 1; // (true)
					break; // break outermost for loop (l)
				}
			}

			if((ok == 0)) {
				printf("ERROR: Missing back link found: Residue name: %s Atom %d:%s --> Atom %d:%s \n",
					atom[j].parentl3name.c_str(),
					atom[j].ir,
					atom[j].rawname.c_str(),
					atom[resti].ir,
					atom[resti].rawname.c_str());
				printf("ERROR: Check Forcefield Definition File \n");
				return_value = -1; // error found
			}
		}
	}
	return return_value;
}

int MoleculeBase::loadParameters(const FFParamSet &ffps)
{
	for(size_t iat=0;iat<nAtoms();iat++)
	{
		int typemol = ffps.findMoleculeType(atom[iat].parentname);
		if(typemol < 0)
		{
			printf("ERROR: Cannot assign Forcefield parameters to atom '%d': '%s', residue '%s':'%d' - residue name unknown\n",
				iat,atom[iat].rawname.c_str(),atom[iat].parentname.c_str(),atom[iat].ir);
			return -1;
		}
		int typeatom = ffps.molecule[typemol].findAtomRaw(atom[iat].rawname);
		if(typeatom < 0)
		{
			printf("ERROR: Cannot assign Forcefield parameters to atom '%d':'%s', residue '%s':'%d' - atom name unknown\n",
				iat,atom[iat].rawname.c_str(),atom[iat].parentname.c_str(),atom[iat].ir);
			return -1;
		}

		atom[iat].CopyParamsFrom(ffps.molecule[typemol].atom[typeatom]);
		atom[iat].FFType = ffps.molecule[typemol].atom[typeatom].FFType; // Atom Type number
	}

	// now that parameters have been loaded detect residue boundaries
	detectSubBoundaries();

	loadedparams = true; // flag that we've successfully executed this function at least once
	return 0;
}

// This function fills the cov13atom and cov14atom arrays
int MoleculeBase::calc1314BondOrders(){
	typedef struct __indextrail{
		int i1, i2, i3, i4; // i1 = index, i2 is i of creator, i3 is creator of creator, etc...
	}indextrail;

	unsigned i, restr;
	int j, listi;
	indextrail list_index[10*4 * 3 * 3 + 1];
	char list_order[10*4 * 3 * 3 + 1];
	int listsize;

	for(i = 0; i < nAtoms(); i++) {
		atom[i].cov13atom.clear();
		atom[i].cov14atom.clear();
	}

	for(i = 0; i < nAtoms(); i++) {

		list_index[0].i1 = i; // add the i atom as first member in list
		list_index[0].i2 = -1; // invalid creator
		list_index[0].i3 = -1; // invalid creator
		list_index[0].i4 = -1; // invalid creator
		list_order[0] = 0; // zero order
		listsize = 1;

		//printf("Processing <%d> ..",i);

		for(listi = 0; listi < listsize; listi++) { // now for every list member
			j = list_index[listi].i1; // get atom index

			//determine its ptr_bondorder and save it *only* if smaller than preexisting
			//value to ensure correctness in cyclic structures ! (atom may be reached from
			// several sides with different bond order - in this case save the lowest!)
			// this restriction is only true for the absolute ptr_bondorder matrix
			// the same end atom can still appear in 1-3 & 1-4 context !
			// hence the above condition only applies to the line below, not
			// to the next block of intructions ( addition to individual atom's lists)

			//also add to individual atoms list (save the index of the atoms)
			switch (list_order[listi]) {
				case 0:
				case 1: // nothing, ince we already have all 1,1 and 1,2 pairs
					break;
				case 2:{
					CovalentAtom newcov13atom;
					newcov13atom.i = list_index[listi].i1; // save 1,3 atom
					newcov13atom.i2 = list_index[listi].i2; // and 1,2 atom
					atom[i].cov13atom.push_back(newcov13atom);
					   }break;
				case 3:{
					CovalentAtom newcov14atom;
					newcov14atom.i = list_index[listi].i1; // save 1,4 atom
					newcov14atom.i2 = list_index[listi].i2; // the 1,3 atom
					newcov14atom.i3 = list_index[listi].i3; // and the 1,2 atom
					atom[i].cov14atom.push_back(newcov14atom);
					   }break;
				default:
					printf("CODE ERROR in calcBondOrders();\n");
					break;
			}

			if(list_order[listi] >= 3)
				continue; // dont add more atoms to list once the 1,4 atom is reached

			for(restr = 0; restr < atom[j].cov12atom.size(); restr++) { // otherwise add further branches
				if(atom[j].cov12atom[restr].i == list_index[listi].i2 )
					continue; //ignore backsteps


				list_index[listsize].i4 = list_index[listi].i3;
				list_index[listsize].i3 = list_index[listi].i2;
				list_index[listsize].i2 = list_index[listi].i1;
				list_index[listsize].i1 = atom[j].cov12atom[restr].i;
				//save atom index of restraint 'restr'
				list_order[listsize] = list_order[listi] + 1; //and increase its ptr_bondorder

				listsize++;
				if(listsize >= (10*4 * 3 * 3)) {
					printf("OVERFLOW ERROR during bond order determination - check code\n");
					return -1;
				}
			}
		}
	}
	return 0;
}

void MoleculeBase::detectSubBoundaries()
{
	atom.trim();
	detectResidueBoundaries();
	detectMoleculeBoundaries();
}

void MoleculeBase::detectResidueBoundaries()
{
	// NOTE - HACK
	// This function currently uses knowledge - which is ~pure~ evil ;-)
	// a more generalised FFPARAM based method is required 

	unsigned ir, i;
	int lastres = -1;
	res.clear();

	// first detect residue boundaries.
	for(size_t i = 0; i < nAtoms(); i++)
	{
		atom[i].i = i;

		if(atom[i].ir != lastres)
		{
			Residue newres;
			newres.ir = atom[i].ir;
			newres.ifirst = i;
			newres.ilast = nAtoms()-1;
			res.push_back(newres);
			lastres = atom[i].ir;
		}	
	}

	// now fill residue special atom indices if possible.
	if(nResidues() > 0 )
	{
		// set last atom (ilast)
		for(ir = 0; (int)ir < int((int)nResidues()-1); ir++) 
		{
			res[ir].ilast = res[ir+1].ifirst - 1;
		}
		res[ir].ilast = nAtoms()-1;

		for(ir = 0; ir < nResidues(); ir++) 
		{
			res[ir].iCA = -1;
			res[ir].iN  = -1;
			res[ir].iC  = -1;
			res[ir].iO  = -1;
			res[ir].iH  = -1;
			res[ir].iCB = -1;
			res[ir].iHA = -1;
		}

		// new & improved loop that takes much less time with large systems
		// find some key atom indices
		for(i = 0;i<nAtoms(); i++){
			if(cmpstring("CA", atom[i].pdbname)){
				if((atom[i].ir>=0) && (atom[i].ir<nResidues())){
					res[atom[i].ir].iCA = i;
				}
			}
			if(cmpstring("N", atom[i].pdbname)){
				if((atom[i].ir>=0) && (atom[i].ir<nResidues())){
					res[atom[i].ir].iN = i;
				}
			}
			if(cmpstring("C", atom[i].pdbname)){
				if((atom[i].ir>=0) && (atom[i].ir<nResidues())){
					res[atom[i].ir].iC = i;
				}
			}
			if(cmpstring("O", atom[i].pdbname)){
				if((atom[i].ir>=0) && (atom[i].ir<nResidues())){
					res[atom[i].ir].iO = i;
				}
			}
			if(cmpstring("H", atom[i].pdbname)){
				if((atom[i].ir>=0) && (atom[i].ir<nResidues())){
					res[atom[i].ir].iH = i;
				}
			}

			if(atom[i].parentletter == 'G') {
				if(cmpstring("HA2", atom[i].pdbname)){
					if((atom[i].ir>=0) && (atom[i].ir<nResidues())){
						res[atom[i].ir].iCB = i;
					}
				}
				if(cmpstring("HA3", atom[i].pdbname)){
					if((atom[i].ir>=0) && (atom[i].ir<nResidues())){
						res[atom[i].ir].iHA = i;
					}
				}
			} else {
				if(cmpstring("CB", atom[i].pdbname)){
					if((atom[i].ir>=0) && (atom[i].ir<nResidues())){
						res[atom[i].ir].iCB = i;
					}
				}
				if(cmpstring("HA", atom[i].pdbname)){
					if((atom[i].ir>=0) && (atom[i].ir<nResidues())){
						res[atom[i].ir].iHA = i;
					}
				}
			}
		}

		for(ir = 0; ir < nResidues(); ir++) 
		{
			res[ir].letter = atom[res[ir].ifirst].parentletter;
			res[ir].param = atom[res[ir].ifirst].parent;
		}
	}
}

void MoleculeBase::detectMoleculeBoundaries()
{
	unsigned imol, i;
	int lastmol = -1;
	mol.clear();

	for(i = 0;i<nAtoms(); i++)
	{
		if(atom[i].imol != lastmol)
		{
			MoleculeRange newmol;
			newmol.name = atom[i].parentname;
			newmol.ifirst = i;
			newmol.ilast = nAtoms()-1;
			mol.push_back(newmol);
			lastmol = atom[i].imol;
		}
	}

	if(nMolecules() > 0 )
	{
		// set last atom (ilast)
		for(imol = 0; (int)imol < int((int)nMolecules()-1); imol++) 
		{
			mol[imol].ilast = mol[imol+1].ifirst - 1;
	
		}
		mol[imol].ilast = nAtoms()-1;
		for(imol = 0; (int)imol < nMolecules() ; imol++) 
		{
			mol[imol].irfirst = atom[mol[imol].ifirst].ir;
			mol[imol].irlast = atom[mol[imol].ilast].ir;
		}
	}
}

void MoleculeBase::append(const MoleculeBase &appendmol)
{
	size_t offset = nAtoms();  // current number of atoms
	size_t res_offset = nResidues(); // current number of residues
	size_t mol_offset = nMolecules(); // current number of (sub)molecules

	for(unsigned iatom=0;iatom<appendmol.nAtoms();iatom++)
	{
		atom.push_back(appendmol.atom[iatom]);
				
		atom[nAtoms()-1].imol += mol_offset;
		atom[nAtoms()-1].ir  += res_offset;

		atom[nAtoms()-1].offsetInternalIndices(offset);
	}
	atom.trim();
	detectSubBoundaries();

	m_Sequence.append( appendmol.m_Sequence );
}

void MoleculeBase::append(const System& appendsys)
{
	// Append molecules	
	size_t res_offset = nResidues(); // current number of residues
	size_t mol_offset = nMolecules(); // current number of (sub)molecules

	for(size_t imol=0;imol<appendsys.nMolecules();imol++)
	{
		size_t offset = nAtoms();  // *current* number of atoms

		const MoleculeBase &appendmol = appendsys.getMolecule(imol);
		for(size_t iatom = 0; iatom < appendmol.nAtoms(); iatom++)
		{
			atom.push_back(appendmol.atom[iatom]);
					
			atom[nAtoms()-1].imol += mol_offset;
			atom[nAtoms()-1].ir += res_offset;
			atom[nAtoms()-1].offsetInternalIndices(offset);
		}

		res_offset += appendmol.nResidues();
		mol_offset += appendmol.nMolecules();

		m_Sequence.append( appendmol.m_Sequence );
	}	
	atom.trim();
	detectSubBoundaries();	
}

void MoleculeBase::beginBonding()
{
	ASSERT( m_BondingMode != true, CodeException, "Bonding mode is already active");
	m_BondingMode = true;
}

void MoleculeBase::endBonding()
{
	ASSERT( m_BondingMode != false, CodeException, "Bonding mode is not already active");
	m_BondingMode = false;
	// Call the calculation of boning orders
	if( 0 != calc1314BondOrders() ) 
		throw( ProcedureException("addBond() caused errors or inconsistencies during the creation of 1:3 and 1:4 covalent structure. \nCheck the forcefield input files and/or input sequences. Ensure that capping \nresidue types are used at the termini.")); 
}

void MoleculeBase::addBond(size_t indexi, size_t indexj )
{
	atom[indexi].add12OnlyCovbond((int)indexj);
	atom[indexj].add12OnlyCovbond((int)indexi);
	// Calculate the overall bonding order if we are not in bonding mode
	if( !m_BondingMode && 0 != calc1314BondOrders() ) 
		throw( ProcedureException("addBond() caused errors or inconsistencies during the creation of 1:3 and 1:4 covalent structure. \nCheck the forcefield input files and/or input sequences. Ensure that capping \nresidue types are used at the termini.")); 
}

void MoleculeBase::addDisulphide(size_t indexi, size_t indexj )
{
	atom[indexi].add12OnlyCovbond((int)indexj);
	atom[indexj].add12OnlyCovbond((int)indexi);
	// Calculate the overall bonding order if we are not in bonding mode
	if( !m_BondingMode && 0 != calc1314BondOrders() ) 
		throw( ProcedureException("addDisulphide() caused errors or inconsistencies during the creation of 1:3 and 1:4 covalent structure. \nCheck the forcefield input files and/or input sequences. Ensure that capping \nresidue types are used at the termini.")); 
}












////////////////////////////////////////////////////////////////////////////////
//
// Coordinate Analysing/Modifying Functions

double MoleculeBase::calcAtomSqrDistance(int i, int j) const
{
	return atomxyz(i).sqrdist(atomxyz(j));
}

double MoleculeBase::calcAtomDistance(int i, int j) const
{
	return atomxyz(i).dist(atomxyz(j));
}

double MoleculeBase::calcAtomAngle(int i, int a, int j) const
{
	dvector b1(atomxyz(i));
	dvector b2(atomxyz(j));
	b1.sub(atomxyz(a));
	b2.sub(atomxyz(a));
	return b1.angleWith(b2);
}

bool MoleculeBase::isMolStartResIndex( int _ResIndex ) const
{
	for( size_t i = 0; i < mol.size(); i++ )
	{
		const MoleculeRange& range = mol.at(i);
		if( range.irfirst == _ResIndex ) return true;
	}
	return false;
}

bool MoleculeBase::isMolEndResIndex( int _ResIndex ) const
{
	for( size_t i = 0; i < mol.size(); i++ )
	{
		const MoleculeRange& range = mol.at(i);
		if( range.irlast == _ResIndex ) return true;
	}
	return false;
}

void MoleculeBase::calcResiduePhiPsi(int ir, double &fphi, double &fpsi) const
{
	dvector i, a, b, j;

	if((ir < nResidues())&&
		(ir > 0)&&
		(res[ir - 1].iC >= 0) &&
		(res[ir].iN >= 0) &&
		(res[ir].iCA >= 0) &&
		(res[ir].iC >= 0) )
	{
		i.setTo(atomxyz(res[ir - 1].iC));
		a.setTo(atomxyz(res[ir].iN));
		b.setTo(atomxyz(res[ir].iCA));
		j.setTo(atomxyz(res[ir].iC));
		fphi = calcTorsionAngle(i, a, b, j);
	} 
	else 
	{
		fphi = 0.0;
	}

	if((ir < (nResidues()-1))&&
		(ir >= 0)&&
		(res[ir].iN >= 0) &&
		(res[ir].iCA >= 0) &&
		(res[ir].iC >= 0) &&
		(res[ir + 1].iN >= 0))
	{
		i.setTo(atomxyz(res[ir].iN));
		a.setTo(atomxyz(res[ir].iCA));
		b.setTo(atomxyz(res[ir].iC));
		j.setTo(atomxyz(res[ir + 1].iN));
		fpsi = calcTorsionAngle(i, a, b, j);
	} 
	else 
	{
		fpsi = 0;
	}
}

void MoleculeBase::calcResidueOmega(int ir, double &fomega)const
{
	dvector i, a, b, j;

	if((ir < nResidues())&&
		(ir > 0)&&
		(res[ir - 1].iCA >= 0) &&
		(res[ir - 1].iC >= 0) &&
		(res[ir].iN >= 0) &&
		(res[ir].iCA >= 0) )
	{
		i.setTo(atomxyz(res[ir - 1].iCA));
		a.setTo(atomxyz(res[ir - 1].iC));
		b.setTo(atomxyz(res[ir].iN));
		j.setTo(atomxyz(res[ir].iCA));
		fomega = calcTorsionAngle(i, a, b, j);
	} 
	else 
	{
		fomega = 0.0;
	}
}

void MoleculeBase::setAllResiduePhiPsi(double nphi, double npsi)
{
	for(unsigned ir=0;ir<nResidues();ir++){
		setResiduePhiPsi(ir,nphi,npsi);
	}
}

void MoleculeBase::setAllResidueOmega(double nomega)
{
	// Begin from 1, as the first residue has no omega angle by definition.
	for(unsigned ir=1;ir<nResidues();ir++){
		setResidueOmega(ir,nomega);
	}
}

void MoleculeBase::setProteinAlphaHelix(){
	setAllResiduePhiPsi(DegToRad(-57.0),DegToRad(-47.0));
	setAllResidueOmega(Maths::MathConst::PI);
}

void MoleculeBase::setProteinReverseHelix(){
	setAllResiduePhiPsi(DegToRad(57.0),DegToRad(47.0));
	setAllResidueOmega(Maths::MathConst::PI);
}

void MoleculeBase::setResiduePhiPsi(int ir, double nphi, double npsi)
{
	double fphi, fpsi;
	ASSERT(((ir>=0)&&(ir<nResidues())), OutOfRangeException, "MoleculeBase::setResiduePhiPsi() out of range");

	calcResiduePhiPsi(ir, fphi, fpsi);

	rotateBond(res[ir].iN, res[ir].iCA, nphi - fphi);
	rotateBond(res[ir].iCA, res[ir].iC, npsi - fpsi);
}

void MoleculeBase::setResidueOmega(int ir, double nomega)
{
	double fomega;
	ASSERT(((ir>0)&&(ir<nResidues())), OutOfRangeException, "MoleculeBase::setResidueOmega() out of range");

	calcResidueOmega(ir, fomega);

	rotateBond(res[ir-1].iC, res[ir].iN, nomega - fomega);
}

int MoleculeBase::setBioSequence( const Sequence::BioSequence& sequenceInformation )
{
	m_Sequence = sequenceInformation; // copy our sequence to the SysMolecule
	return 0;
}

const Sequence::BioSequence& MoleculeBase::getSequence() const
{
	return m_Sequence;
}

std::string MoleculeBase::getResName( size_t _ir ) const
{
	return m_Sequence.getResidue(_ir).getResName();
}

std::string MoleculeBase::getFullResName( size_t _ir ) const
{
	return m_Sequence.getResidue(_ir).getResName();
}

void MoleculeBase::resetAllAtomMovedFlags()
{
  for(int i = 0; i < nAtoms(); i++) {
    atom[i].setMoved(false);
  }
}

double MoleculeBase::getTotalMass() const
{
	double totalmass = 0;
	for(int i = 0; i < nAtoms(); i++) {
		totalmass += atom[i].mass;
	}
	return totalmass;
}

double MoleculeBase::getTotalCharge() const
{
	double totalcharge = 0;
	for(int i = 0; i < nAtoms(); i++) {
		totalcharge += atom[i].charge;
	}
	return totalcharge;
}


dvector MoleculeBase::getCentreOfMass() const
{
	dvector com(0,0,0);
	dvector p;
	for(int i = 0; i < nAtoms(); i++) {
		p.setTo(atomxyz(i));
		p.mul(atom[i].mass);
		com.add(p);
	}
	com.div(getTotalMass());
	return com;
}

dvector MoleculeBase::getCentreOfAvailableGeometry() const
{
	dvector cog(0,0,0);
	int cnt = 0;
	for(int i = 0; i < nAtoms(); i++)
	{
		if( 
			//atom[i].isValid() &&
			!atom[i].isRebuildRequired() &&			
			!atom[i].isDummy() 
			)
		{
			cnt++;
			cog.add(atomxyz(i));
		}
	}
	if( cnt != 0 ) cog.div((double)cnt);
	return cog;
}

dvector MoleculeBase::getCentreOfGeometry() const
{
	dvector cog(0,0,0);
	for(int i = 0; i < nAtoms(); i++)
		cog.add(atomxyz(i));
	cog.div((double)nAtoms());
	return cog;
}

dvector MoleculeBase::getCentreOfGeometry(unsigned i) const 
{
	return getCentreOfGeometry(mol[i].ifirst,mol[i].ilast);
}

dvector MoleculeBase::getCentreOfGeometry(const PickBase& _picker) const
{
	dvector cog(0.0,0.0,0.0);
	double count = 0.0;
	for(int i = 0; i < nAtoms(); i++)
		if( _picker.matches(atom[i]) )
		{
			cog.add(atomxyz(i));
			count+=1.0;
		}
	if( count != 0.0 )
		cog.div(count);
	return cog;
}

dvector MoleculeBase::getCentreOfGeometry(unsigned ifirst, unsigned ilast) const 
{
	dvector cog(0,0,0);
	for(int i = ifirst; i <= ilast; i++){
		cog.add(atomxyz(i));
	}
	cog.div((double)(ilast-ifirst+1));
	return cog;
}

void MoleculeBase::calcInertiaTensor(matrix3x3 & I) const{
	matrix3x3 Ii;
	I.setToNull();
	dvector posSI;
	dvector com = getCentreOfMass();
	for(int i = 0; i < nAtoms(); i++) {
		posSI.setTo(atomxyz(i));
		posSI.sub(com);
		posSI.mul(Physics::PhysicsConst::Angstrom);
		Ii.setTo( sqr(posSI.y) + sqr(posSI.z), -posSI.x * posSI.y, -posSI.x * posSI.z,
			-posSI.y * posSI.x ,  sqr(posSI.x) + sqr(posSI.z),	-posSI.y * posSI.z,
			-posSI.z * posSI.x , -posSI.z * posSI.y, sqr(posSI.x) + sqr(posSI.y));
		Ii.mul(atom[i].mass);
		I.add(Ii);
	}
}

void MoleculeBase::printInertiaInfo() const{
	matrix3x3 I;
  calcInertiaTensor(I);
	printf("Ixx =        %.6e   kg.m^2 \n", I.r[0][0] ); 
	printf("Iyy =        %.6e   kg.m^2 \n", I.r[1][1] ); 
	printf("Izz =        %.6e   kg.m^2 \n", I.r[2][2] ); 
	printf("Ixy = Iyx =  %.6e   kg.m^2 \n", I.r[0][1] ); 
	printf("Ixz = Izx =  %.6e   kg.m^2 \n", I.r[0][2] ); 
	printf("Izy = Iyz =  %.6e   kg.m^2 \n", I.r[2][1] ); 
  printf("I = IA IB IC =  %.6e  (kg.m^2)^3 \n",
           I.r[0][0]*I.r[1][1]*I.r[2][2] +
           2.0 * I.r[0][1]*I.r[0][2]*I.r[2][1] -
           I.r[0][0] * sqr(I.r[2][1]) -
           I.r[1][1] * sqr(I.r[2][0]) -
           I.r[2][2] * sqr(I.r[0][1]) );
}

double MoleculeBase::calcRotationalPartition(double temp, unsigned SymNumber) const{
  matrix3x3 I;
  calcInertiaTensor(I);

  // product of principal moments of inertia
  double IAIBIC =
           I.r[0][0]*I.r[1][1]*I.r[2][2] +
           2.0 * I.r[0][1]*I.r[0][2]*I.r[2][1] -
           I.r[0][0] * sqr(I.r[2][1]) -
           I.r[1][1] * sqr(I.r[2][0]) -
           I.r[2][2] * sqr(I.r[0][1]);
	double qrot;
	qrot = 8.0*sqr(MathConst::PI)/(double(SymNumber));
	qrot *= sqrt(cube(2.0*MathConst::PI * Physics::PhysicsConst::kB * temp/sqr(Physics::PhysicsConst::planck))*IAIBIC);
	return qrot;
}

void MoleculeBase::setToGeometry(){

	for(size_t i = 0; i < nAtoms(); i++) 
	{
		atomxyz(i).setTo(atom[i].posGeom());
	}
}

void MoleculeBase::setToReference(){

	for(size_t i = 0; i < nAtoms(); i++) 
	{
		atomxyz(i).setTo(atom[i].posRef());
	}
}

void MoleculeBase::setAsReference(){

	for(size_t i = 0; i < nAtoms(); i++) 
	{
		atom[i].posRef().setTo(atom[i].pos());
	}
}


double MoleculeBase::calcCRMS_AllAtom( bool useOnlyKnownAtom ) const{
	double cRMS;
	int i;
	int hAtoms = 0;

	dvector *curatompos = new dvector[nAtoms()];
	dvector *tgtatompos = new dvector[nAtoms()];

	if( useOnlyKnownAtom )
	{
		for(i = 0; i < nAtoms(); i++) 
		{
			if( atom[i].isKnownStructure() )
			{
				curatompos[hAtoms].setTo(atomxyz(i));
				tgtatompos[hAtoms].setTo(atom[i].posRef());
				hAtoms++;
			}
		}
	}
	else
	{
		hAtoms = (int) nAtoms();
		for(i = 0; i < nAtoms(); i++) 
		{
			curatompos[hAtoms].setTo(atomxyz(i));
			tgtatompos[hAtoms].setTo(atom[i].posRef());
		}
	}

	cRMS = calcVectorCRMS(curatompos, tgtatompos, hAtoms);
	delete[]curatompos;
	delete[]tgtatompos;
	return cRMS;
}

double MoleculeBase::calcCRMS_HeavyAtom( bool useOnlyKnownAtom ) const{
	double hacRMS;
	int i;
	int hAtoms = 0;

	dvector *curatompos = new dvector[nAtoms()];
	dvector *tgtatompos = new dvector[nAtoms()];

	if( useOnlyKnownAtom )
	{
		for(i = 0; i < nAtoms(); i++) 
		{
			if(atom[i].isHydrogen() || !atom[i].isKnownStructure()) continue; // exclude hydrogens and unknown structure
			curatompos[hAtoms].setTo(atomxyz(i));
			tgtatompos[hAtoms].setTo(atom[i].posRef());
			hAtoms++;
		}
	}
	else
	{
		for(i = 0; i < nAtoms(); i++) 
		{
			if(atom[i].isHydrogen()) continue; // exclude hydrogens
			curatompos[hAtoms].setTo(atomxyz(i));
			tgtatompos[hAtoms].setTo(atom[i].posRef());
			hAtoms++;
		}
	}

	hacRMS = calcVectorCRMS(curatompos, tgtatompos, hAtoms);
	delete[]curatompos;
	delete[]tgtatompos;
	return hacRMS;
}

double MoleculeBase::calcCRMS_CA( bool useOnlyKnownAtom ) const
{
	double CacRMS;
	int hAtoms = 0;

	dvector *curatompos = new dvector[nAtoms()];
	dvector *tgtatompos = new dvector[nAtoms()];

	if( useOnlyKnownAtom )
	{
		for(size_t i = 0; i < nAtoms(); i++) 
		{
			if(!atom[i].isCAlpha() || !atom[i].isKnownStructure() ) continue; // only include CA atoms
			curatompos[hAtoms].setTo(atomxyz(i));
			tgtatompos[hAtoms].setTo(atom[i].posRef());
			hAtoms++;
		}
	}
	else
	{
		for(size_t i = 0; i < nAtoms(); i++) 
		{
			if(!atom[i].isCAlpha()) continue; // only include CA atoms
			curatompos[hAtoms].setTo(atomxyz(i));
			tgtatompos[hAtoms].setTo(atom[i].posRef());
			hAtoms++;
		}
	}

	CacRMS = calcVectorCRMS(curatompos, tgtatompos, hAtoms);

	delete[]curatompos;
	delete[]tgtatompos;

	return CacRMS;
}

void MoleculeBase::zeroCentreOfAvailableGeometry()
{
	dvector cog = getCentreOfAvailableGeometry();
	for(size_t i = 0; i < nAtoms(); i++) 
	{
		atomxyz(i).sub(cog);
	}
}

void MoleculeBase::zeroCentreOfGeometry()
{
	dvector cog = getCentreOfGeometry();
	for(size_t i = 0; i < nAtoms(); i++) 
	{
		atomxyz(i).sub(cog);
	}
}

void MoleculeBase::zeroCentreOfMass()
{
	dvector com = getCentreOfMass();
	for(size_t i = 0; i < nAtoms(); i++) 
	{
		atomxyz(i).sub(com);
	}
}


void MoleculeBase::alignAlongPrincipalAxes() 
{
	matrix3x3 I;
	calcInertiaTensor(I);
	double lambda1,lambda2,lambda3;
	dvector p1,p2,p3;
	I.diagonaliseSymetric(lambda1,lambda2,lambda3,p1,p2,p3);
	matrix3x3 rot;
	dvector i(1,0,0);
  dvector j(0,1,0);
	// rotate the system so that the first tow principal 
	// components point along the cartesian axe
  // by definition the last will also point inthe right direction
	superimpose(i,j,p1,p2,rot);   
	rotate(rot);
}


void MoleculeBase::moveMolecule(int imol, const Maths::dvector & disp)
{
	moveParticles(mol[imol].ifirst, mol[imol].ilast+1, disp);
}

void MoleculeBase::moveParticles(size_t start, size_t end, const dvector & disp)
{
	for(size_t i = start; i < end; i++)
		atomxyz(i).add(disp);
}

void MoleculeBase::rotateMolecule(int imol, const Maths::dvector & centre, const Maths::matrix3x3 & rmat)
{
	rotateParticles(mol[imol].ifirst, mol[imol].ilast, centre, rmat);
}

void MoleculeBase::rotateParticles(size_t start, size_t end, const dvector & centre, const matrix3x3 & rmat)
{
	for(size_t i = start; i < end; i++) 
	{
		atomxyz(i).sub(centre);
		atomxyz(i).mulmat(rmat);
		atomxyz(i).add(centre);
	}
}


int MoleculeBase::rotateBond(int bi, int bj, double angle){
	// first check atom indices;
	if(bi<0) return -1;
	if(bj<0) return -1;

	int n;
	matrix3x3 rmat;
	dvector axis;

	// first check that this is actually a bond - by finding bj on bi's covalent list
	bool foundBJ = false;
	for(int i = 0; i < atom[bi].cov12atom.size(); i ++ ){
		if(atom[bi].cov12atom[i].i == bj)
		{
			foundBJ = true;
			break;
		}
	}
	if(!foundBJ) return -1; /// if we didnt find it and the for loop ran out - abort

	axis.setTo(atomxyz(bj));
	axis.sub(atomxyz(bi));

	rmat.setToAxisRot(axis, angle);

	int* localStatus = new int[nAtoms()];
	for(int i = 0; i < nAtoms(); i++) {
		localStatus[i] = -1; // mark as unrotated
	}

	int *list = new int[nAtoms()];
	int listentries = 1;
	int curentry = -1;
	int listatom;
	int neighatom;

	localStatus[bi] = 2; // mark i as rotated (so we cant traverse the i-j bond
	localStatus[bj] = 1; // mark j as listed
	list[0] = bj; // set j as the first entry

	while(curentry < (listentries - 1)) {
		curentry++;

		// first rotate the list member:

		listatom = list[curentry];

		if(localStatus[listatom] >= 2)
			continue; // ignore if already rotated

		atomxyz(listatom).sub(atomxyz(bi)); // translate to axis origin
		atomxyz(listatom).mulmat(rmat); // rotate
		atomxyz(listatom).add(atomxyz(bi)); // translate back to original position
		localStatus[listatom] = 2; // mark as rotated

		// now find listatom's covalently linked atoms to create new entries to the list

		for(n = 0; n < atom[listatom].cov12atom.size(); n++) {
			neighatom = atom[listatom].cov12atom[n].i; // get atom index of cov neighbour
			if(neighatom == listatom)
				continue; // ignore direct backstepping
			if(localStatus[neighatom] >= 1)
				continue; // ignore if already listed

			localStatus[neighatom] = 1; // mark as listed

			// if arive here - add to list for processing
			list[listentries] = neighatom;
			listentries++;
			if(listentries >= nAtoms()) {
				THROW(CodeException,"FATAL ERROR: List overflow in MoleculeBase::rotateBond(...) \n");
			}
		}
	}

	delete[]list; // free list memory
	delete[]localStatus;

	return 0;
}



// gets maximum absolute coordinate in system (i.e. largest value of |x|,|y| and |z|
Maths::dvector MoleculeBase::getEncompassingVector() const
{
	Maths::dvector largest(0,0,0);
	int i;
	for(i = 0; i < nAtoms(); i++) {
		if(fabs(atomxyz(i).x) > largest.x ) largest.x = fabs(atomxyz(i).x);
		if(fabs(atomxyz(i).y) > largest.y ) largest.y = fabs(atomxyz(i).y);
		if(fabs(atomxyz(i).z) > largest.z ) largest.z = fabs(atomxyz(i).z);
	}
	return largest;
}

// gets maximum coordinate in system (i.e. largest value of |x|,|y| and |z|
Maths::dvector MoleculeBase::getMaximumVector() const
{
	Maths::dvector largest(-DBL_MAX,-DBL_MAX,-DBL_MAX);
	int i;
	for(i = 0; i < nAtoms(); i++) {
		if(atomxyz(i).x > largest.x ) largest.x = atomxyz(i).x;
		if(atomxyz(i).y > largest.y ) largest.y = atomxyz(i).y;
		if(atomxyz(i).z > largest.z ) largest.z = atomxyz(i).z;
	}
	return largest;
}


// gets minimum coordinate in system (i.e. smallest value of |x|,|y| and |z|
Maths::dvector MoleculeBase::getMinimumVector() const
{
	Maths::dvector smallest(DBL_MAX,DBL_MAX,DBL_MAX);
	int i;
	for(i = 0; i < nAtoms(); i++) {
		if(atomxyz(i).x < smallest.x ) smallest.x = atomxyz(i).x;
		if(atomxyz(i).y < smallest.y ) smallest.y = atomxyz(i).y;
		if(atomxyz(i).z < smallest.z ) smallest.z = atomxyz(i).z;
	}
	return smallest;
}

void MoleculeBase::printPDB(const std::string& targetFile)
{
	IO::PDB_Writer rw( targetFile, false );
	rw.write(*this);
}

void MoleculeBase::printPDB()
{
	IO::PDB_Writer rw( false );
	rw.write(*this);
}

size_t MoleculeBase::memuse(int level)
{
	level--;
	size_t s_self=sizeof(MoleculeBase);
	size_t s_atom = 0;
	for(size_t i = 0; i< nAtoms(); i++) s_atom += atom[i].memuse(level);
	size_t s_res = 0;
	for(size_t i = 0; i< nResidues(); i++) s_res += res[i].memuse(level);
	size_t s_mol = 0;
	for(size_t i = 0; i< nMolecules(); i++) s_mol += mol[i].memuse(level);
	if(level >= 0){
		for(int i=0;i<level;i++) printf(" ");
		printf("MoleculeBase:  %d : %d %d %d %d\n",
		                    (int)(s_self+ s_atom+ s_res+ s_mol),
				(int)s_self, (int)s_atom, (int)s_res, (int)s_mol);
	}
	return s_mol + s_res + s_atom + s_self;
}



bool operator==(const MoleculeBase &mol1, const MoleculeBase &mol2 ){
	if(mol1.nAtoms()!=mol2.nAtoms()) return false;
	if(mol1.nResidues()!=mol2.nResidues()) return false;
	if(mol1.nMolecules()!=mol2.nMolecules()) return false;
	if(!(mol1.m_Sequence==mol2.m_Sequence)) return false;
	for(int i = 0; i < mol1.nAtoms(); i++) {
		if(mol1.atom[i].rawname != mol2.atom[i].rawname) return false;
		if(mol1.atom[i].pdbname != mol2.atom[i].pdbname) return false;
	}
	return true; 
} 

void Molecule::save( IO::OutputFile &_output ){
	_output.save( *this );
}

