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

// mmLib includes
#include "sequence/sequence.h"
#include "system/molecule.h"

#include "system/genpolymer.h"

using namespace Maths;
using namespace Library;

int GeneratorCore( const Sequence::BioSequence &bioSeq, Molecule &newmol, const FFParamSet &ffps, bool _Polymerise, Verbosity::Type OutputLevel)
{
	int nres = bioSeq.size(); // number of residues
	int i, ir, iat, nnewparts;
	std::string iresname;

	newmol.name = "Polymer";

	std::vector<const MoleculeDefinition *> resdef; // array of pointers to the forcefield defs

	for(i = 0; i < nres; i++) {
		iresname = bioSeq.getFullName(i);
		if(iresname.compare("???")==0)
		{
			printf("ERROR: Unknown residue passed in bioseq to GeneratePolymer()!\n");
			return -1;
		}
		// Sucessfully isolated identifier, now find the molecule
		int imolecule = -1;
		imolecule = ffps.findMoleculeType(iresname);
		if(imolecule < 0) {
			printf("ERROR: Residue name/alias '%s' unknown! \n", iresname.c_str());
			return -1;
		}
		resdef.push_back(&ffps.molecule[imolecule]);
	}

	// Work out how many particles are needed
	nnewparts = 0;
	for(ir = 0; ir < resdef.size(); ir++)
		nnewparts += resdef[ir]->atom.size();

	if(OutputLevel)
		printf("Adding %d particles, %d residues to system ...\n", nnewparts, nres);

	for(ir = 0; ir < resdef.size(); ir++)
    {
		if(OutputLevel)
        {
			if((ir>0)&&((ir%15)==0)) printf("\n");
			printf("%4s ", resdef[ir]->c_name());
		}

		if(ir == 0)
        { 
            // the first residue must have no back links
            if(resdef[ir]->hasBackLink())
            {
				printf("WARNING: The first residue is not a start-of-polymer residue \n");
			}
		}
		if(ir == (resdef.size() - 1))
        { 
            // the last residue must have no forward links
            if(resdef[ir]->hasFrwdLink())
            {
				printf("WARNING: The last residue is not a end-of-polymer residue \n");
			}
		}

		for(iat = 0; iat < resdef[ir]->atom.size(); iat++)
        {
			Particle newatom(resdef[ir]->atom[iat]); // create a 'copy' of the template atom
			newatom.ir = ir; // residue number
			newatom.parentname = newatom.parent->s_name();
			newatom.parentl3name = newatom.parent->s_l3_name();
			newatom.parentletter = newatom.parent->letter;
			newmol.addParticle(newatom);
		}
	}
	if(OutputLevel) printf("\n\n");

	// create residue container by analysing the atom array
	newmol.detectSubBoundaries();
	if(OutputLevel) printf("Detected residue boundaries\n");

	if( _Polymerise )
	{
		// now move each residue into it's respective position
		if(OutputLevel) printf("Performing Polymerisation: ");
		if(newmol.res.size() != resdef.size())
        {
			printf("\nCODE ERROR: Numbers of residues does not agree\n");
			return -1;
		}

		// Polymerise our new atomic structure from the 'resdef' forward and back link definititions
		for(int ir = resdef.size()-1; ir > 0; ir--)
		{
			int iBackAtom = newmol.findParticle(ir-1,resdef[ir]->backName);
			int iFrwdAtom = newmol.findParticle(ir,resdef[ir-1]->frwdName);

			if(iBackAtom < 0)
            {
				printf("\nERROR: Cannot find BackLink Atom in previous residue of residue %d\n", ir);
				return -1;
			}
			if(iFrwdAtom < 0)
            {
				printf("\nERROR: Cannot find Forward Link Atom in next residue of residue %d\n", ir);
				return -1;
			}

			dvector resoffset;
			matrix3x3 resrot;
			superimpose(
				newmol.atom[iBackAtom].pos(),
				resdef[ir-1]->frwdpos,
				resdef[ir]->backpos,
				newmol.atom[iFrwdAtom].pos(),
				resrot);

			dvector subt(resdef[ir]->backpos);
			dvector addt(newmol.atom[iBackAtom].pos());

			for(size_t i = newmol.res[ir].ifirst; i < newmol.nAtoms(); i++)
			{
				D_ASSERT(newmol.atom[i].ir >= ir,CodeException,"Index logic error!");
				newmol.atom[i].pos().sub(subt);
				newmol.atom[i].pos().mulmat(resrot);
				newmol.atom[i].pos().add(addt);
			}

			dvector backpos;
			backpos = resdef[ir]->backpos;
			backpos.sub(subt);
			backpos.mulmat(resrot);
			backpos.add(addt);
		}

		if(OutputLevel) printf(" Done!\n");
	}

	if(OutputLevel) printf("Covalent structure:\n  Creating: ");
	newmol.createCovalentStructure();
	// check if all bonds have backbonds, i.e. check integrity of covalent structure
	if(OutputLevel) printf("Done!\n  Checking: ");
	if(newmol.checkCovalentStructure() != 0) return -1;
	if(OutputLevel) printf("Done!\n");

	// Try this, doesnt matter if it doesnt work - it will fail silently, and it doesn't matter if it does,
	// but it will give us a nicer helical representation as the geometric / default reference state.
	newmol.setProteinAlphaHelix();

	// Set atom flags and geometrical positions
	for( size_t i = 0; i < newmol.atom.size(); i++ )
	{
		newmol.atom[i].setKnownStructure(false); // this structure is not from atomic data
		newmol.atom[i].setRebuildRequired(false); // No building is required - the structure should be sensible (its made from 'idealised' FF definitions). Other routines can override this later ...
		newmol.atom[i].posGeom().setTo(newmol.atom[i].pos()); // Assign our Geometric structure is current reference structure
	}

	// target structure is current reference structure (this is overridden by file imports)
	newmol.setAsReference();

	if(OutputLevel) printf("Loading MoleculeBase Parameters\n");
	newmol.loadParameters(ffps);

	// Assign the sequence we have just built
	newmol.setBioSequence(bioSeq);

	if(OutputLevel) printf("Complete!\n");
	return 0;
}

int GeneratePolymer( const std::string &aseq, Molecule &newmol, const FFParamSet &ffps, Verbosity::Type OutputLevel)
{
	Sequence::BioSequence bioSeq(ffps); // Create a new sequence using the ffps as the name-mapper
	bioSeq.setNamingException(true);
	bioSeq.append(aseq);
	return GeneratorCore( bioSeq, newmol, ffps, true, OutputLevel );
}

int GeneratePolymer( const Sequence::BioSequence &bioSeq, Molecule &newmol, const FFParamSet &ffps, Verbosity::Type OutputLevel)
{
	return GeneratorCore( bioSeq, newmol, ffps, true, OutputLevel );
}

int GenerateMoleculeSet( const std::string &aseq, Molecule &newmol, const FFParamSet &ffps, Verbosity::Type OutputLevel)
{
	Sequence::BioSequence bioSeq(ffps); // Create a new sequence using the ffps as the name-mapper
	bioSeq.setNamingException(true);
	bioSeq.append(aseq);
	return GeneratorCore( bioSeq, newmol, ffps, false, OutputLevel );
}

int GenerateMoleculeSet( const Sequence::BioSequence &bioSeq, Molecule &newmol, const FFParamSet &ffps, Verbosity::Type OutputLevel)
{
	return GeneratorCore( bioSeq, newmol, ffps, false, OutputLevel );
}

int loadMolecule(const Sequence::ResidueInfo& _Res, Molecule &newmol, const FFParamSet &ffps)
{
	int result = loadMolecule( _Res.getFullName(), newmol,ffps);
	if( result != 0 ) return result;
	Sequence::BioSequence seq;
	seq.append( _Res );
	newmol.setBioSequence( seq );
	return 0;
}

int loadMolecule(const std::string &MolName, Molecule &newmol, const FFParamSet &ffps)
{
	unsigned i, k;
	int molnumber;

	molnumber = ffps.findMoleculeType(MolName);
	if(molnumber < 0) {
		printf("ERROR: Molecule Type '%s' unknown \n", MolName.c_str());
		return -1;
	}
	if(ffps.molecule[molnumber].atom.size() <= 0) {
		printf("ERROR: Molecule '%s' has 0 atoms - cannot add molecule, check definition files\n", MolName.c_str());
		return -1;
	}

	newmol = Molecule(ffps); // clear the molecule container
	newmol.name = MolName;

	for(i = 0; i < ffps.molecule[molnumber].atom.size(); i++)
	{
		Particle newatom(ffps.molecule[molnumber].atom[i]);

		//newatom.ir = ir; // residue number
		newatom.parentname = newatom.parent->s_name();
		newatom.parentl3name = newatom.parent->s_l3_name();
		newatom.parentletter = newatom.parent->letter;
		newmol.addParticle(newatom);
	}

	newmol.beginBonding();
	for(i = 0; i < ffps.molecule[molnumber].atom.size(); i++)
	{
		// add atom-restraints to restraint list
		for(k = 0; k < ffps.molecule[molnumber].atom[i].r_covalent.size(); k++)
		{
			// go through all covalent restraints
			newmol.addBond( i, ffps.molecule[molnumber].atom[i].r_covalent[k].i );
		}
	}
	newmol.endBonding();

	// Create Residue Array - for a predefined molecule this will merely give 1 residue of course
	newmol.detectSubBoundaries();
	newmol.setAsReference(); // target structure is starting structure
	if( 0 != newmol.calc1314BondOrders() ) return -1;

	newmol.m_Sequence = Sequence::BioSequence(ffps);
	newmol.m_Sequence.append(MolName);

	return 0;
}

// Public interface functions

Molecule NewProtein(FFParamSet &ffps, const std::string &aseq)
{
	Molecule newmol(ffps);
	Molecule emptymol(ffps);
	if(GeneratePolymer(aseq,newmol,ffps,Verbosity::Silent)!=0)
	{
		throw(ProcedureException("Error occured during protein generation"));
	}
	if(newmol.nResidues() < 20) newmol.name = "Peptide";
	else                   newmol.name = "Protein";
	newmol.setAllResiduePhiPsi(Maths::MathConst::PI,Maths::MathConst::PI);
	newmol.setAllResidueOmega(Maths::MathConst::PI);
	newmol.setAsReference();
	return newmol;
}


Molecule NewProteinHelix(FFParamSet &ffps,const std::string &aseq){
	Molecule newmol(ffps);
	if(GeneratePolymer(aseq,newmol,ffps,Verbosity::Normal)!=0){
		throw(ProcedureException("Error occured during protein generation"));
	}
	if(newmol.nResidues() < 20) newmol.name = "Peptide";
	else                   newmol.name = "Protein";
	newmol.setAllResiduePhiPsi(DegToRad(-57.0),DegToRad(-47.0));
	newmol.setAllResidueOmega(Maths::MathConst::PI);
	newmol.setAsReference();
	return newmol;
}

Molecule NewMolecule(FFParamSet &ffps, const std::string &MolName)
{
	Molecule newmol(ffps);
	if(loadMolecule(MolName,newmol,ffps)!=0){
		Molecule emptymol(ffps);
		fprintf(stderr,"ERROR: Error during loading of molecule '%s'\n",MolName.c_str());
		throw(ProcedureException("Error occured during molecule generation"));
	}
	newmol.setAsReference();
	return newmol;
}

