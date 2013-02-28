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

#include "forcefields/ffbonded.h"

#include "fileio/pdb.h"
#include "workspace/workspace.h"

// namespace includes
using namespace Maths;

namespace Physics{

	void FF_Bonded::setup()
	{ 
		if(OutputLevel) printf("FF_Bonded Setup:\n\tBonds... \n");
		if(assembleBondList() != 0) throw(ProcedureException("Error occured during setup of Bond List")); 
		if(OutputLevel) printf("\tAngles... \n");
		if(assembleAngleList() != 0) throw(ProcedureException("Error occured during setup of Angle List")); 
		if(OutputLevel) printf("\tTorsions... \n");
		if(assembleTorsionList() != 0) throw(ProcedureException("Error occured during setup of Torsion List")); 
		if(OutputLevel) printf("\tImpropers... \n");
		if(assembleImproperList() != 0) throw(ProcedureException("Error occured during setup of Improper List")); 

		needsetup = false;
	}

	/////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////
	/// \brief Assembles a list of all 1,2 bonds in the system by collecting all 1,2 links from
	/// all atoms and recording in bond[] and bond.size()
	///
	/// Called internally by setup()
	///
	/// \author M.Tyka

	int FF_Bonded::assembleBondList()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		ParticleStore& atomparam = wspace.atom;
		SnapShotAtom *atom = wspace.cur.atom;

		int c, bindex;
		size_t totalbonds;
		int guesscount = 0;

		if(!wspace.ffps().autobonds)
			return 0;

		bond.clear(); //set counters to zero
		totalbonds = 0;

		for( size_t i = 0; i < wspace.atom.size(); i++) //count total number of bonds
			totalbonds += atomparam[i].cov12atom.size();

		if((totalbonds % 2) != 0) {
			THROW(CodeException,"ERROR: Total number of covalent 1,2 links not even - suspect code error/def. file error!");
		}
		totalbonds /= 2; //half total as every bond is only recorded once!

		for(size_t i = 0; i < wspace.atom.size(); i++) {
			for(c = 0; c < atomparam[i].cov12atom.size(); c++) { //collect all the covalent links, assume they're all bond...

				Bond newbond;

				if(i < atomparam[i].cov12atom[c].i)
				{
					// only record half (i.e. don't doubly record 1-3 and 3-1)
					continue; 
				}

				newbond.i = (int)i; //record atoms
				newbond.j = atomparam[i].cov12atom[c].i;

				if(!IgnoreForcefieldParams){
					//find bond Type
					bindex = wspace.ffps().findBondType(atomparam[newbond.i].FFType, 
																			atomparam[newbond.j].FFType);

					//printf(" %s %s \n",);

					if(bindex >= 0) { // An angletype was found successfully, load parameters
						newbond.k = wspace.ffps().BondType[bindex].forceconstant; // force constant
						newbond.l = wspace.ffps().BondType[bindex].length; // eq. bond length

					} else { // otherwise deduce equilibrium length from inital structure, & give guessed k
						newbond.k = 360 * PhysicsConst::kcal2J / PhysicsConst::Na; //force constant
						newbond.l = wspace.calcAtomDistance(newbond.i, newbond.j); //bond length
						printf("WARNING: Bond Type %s-%s not found, guessing parameters (k=720kcal/mol, d=%f.2)\n", 
							wspace.ffps().AtomType[atomparam[newbond.i].FFType].name.c_str(),
							wspace.ffps().AtomType[atomparam[newbond.j].FFType].name.c_str(),
							newbond.l);
						guesscount ++;

						if(guesscount > 10){
							throw( ProcedureException("Too many bond types were not found (see warnings above)\n\
check the forcefield definition files and add the missing parameters!") );
						}
					}
				}else{
						// just guess anyway (user asked us to ignore the forcefield parameters)
						newbond.k = 360 * PhysicsConst::kcal2J / PhysicsConst::Na; //force constant
						newbond.l = wspace.calcAtomDistance(newbond.i, newbond.j); //bond length
				}

				bond.push_back(newbond);
			}
		}

		if(totalbonds != bond.size()) 
		{
			THROW(CodeException,"CODE ERROR #1 in assembleBondList()");
		}

		for(size_t i = 0; i < wspace.atom.size(); i++) 
		{ 
			// now re-collect the equilibrum distances into cov12atom
			for(c = 0; c < atomparam[i].cov12atom.size(); c++) 
			{
				for(bindex = 0; bindex < bond.size(); bindex++) 
				{ 
					// find that bond in bondlist
					if(((bond[bindex].i == i) &&
						(bond[bindex].j == atomparam[i].cov12atom[c].i)) ||
						((bond[bindex].j == i) && (bond[bindex].i == atomparam[i].cov12atom[c].i)))
						break; // finish for loop !
				}
				if(bindex >= bond.size()) 
				{
					THROW(CodeException,"CODE ERROR #2 in assembleBondList()");
				}
				atomparam[i].cov12atom[c].d = bond[bindex].l; // record equilibrium length :)
			}
		}
		return 0;
	}




	// assembleAngleList(); ---------------------------------------------------------
	//
	// Assembles a list of all 1,3 angles in the system by collecting all 1,3 links from
	// all atoms and recording in angle[] and angle.size()
	// ---------------------------------------------------------------------------------

	int FF_Bonded::assembleAngleList()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		ParticleStore& atomparam = wspace.atom;
		SnapShotAtom *atom = wspace.cur.atom;

		int aindex;

		if(!wspace.ffps().autoangles)
			return 0;

		angle.clear(); //set counters to zero
		size_t totalangles = 0;

		for(size_t i = 0; i < wspace.atom.size(); i++) // count total number of angless
		{
			totalangles += atomparam[i].cov13atom.size();
		}

		if((totalangles % 2) != 0) {
			printf("ERROR: Total number of covalent 1,3 links not even - suspect code error/def. file error!\n");
			return -1;
		}
		totalangles /= 2; //half total as every angle is only recorded once!

		//printf("Angle number: %d \n",totalangles);
		for(size_t i = 0; i < wspace.atom.size(); i++) {
			//printf(" Atom %d %d\n",i,atomparam[i].n_cov13atoms);

			for(size_t c = 0; c < atomparam[i].cov13atom.size(); c++) { //collect all the covalent links, assume they're all angle...
				if(i < atomparam[i].cov13atom[c].i)
					continue;//only record half (i.e. don't double record 1-3 and 3-1)
				//printf("atom %d\n",i);

				Angle newangle;
				newangle.i = (int) i; //record atoms
				newangle.a = atomparam[i].cov13atom[c].i2;
				newangle.j = atomparam[i].cov13atom[c].i;

				if(newangle.i < 0){
					THROW(CodeException,"FF_Bonded::assembleAngleList() failure (0)!");
				}
				if(newangle.a < 0){
					THROW(CodeException,"FF_Bonded::assembleAngleList() failure (1)!");
				}
				if(newangle.j < 0){
					THROW(CodeException,"FF_Bonded::assembleAngleList() failure (2)!");
				}

				if(!IgnoreForcefieldParams){

					//find angle Type
					aindex = wspace.ffps().findAngleType(atomparam[newangle.i].FFType,
						atomparam[newangle.a].FFType, atomparam[newangle.j].FFType);
					if(aindex >= 0) { // An angletype was found successfully, load parameters
						newangle.k = wspace.ffps().AngleType[aindex].forceconstant; // those two are defined by the
						newangle.theta0 = wspace.ffps().AngleType[aindex].angle; // ff parameter file

						newangle.l = // this one isnt so far
							wspace.calcAtomDistance(newangle.i, newangle.j); // angle length


						//printf("P %d %lf \n",angle.size(),newangle.k * PhysicsConst::J2kcal * PhysicsConst::Na);
					} else {
						printf("WARNING: In aa %d(atoms %d %d %d) Angle Type not found, guessing parameters \n",
							atomparam[newangle.i].ir, newangle.i, newangle.a, newangle.j);

						newangle.k = 0.03 * PhysicsConst::kcal2J / PhysicsConst::Na; //force constant per deg

						newangle.l = wspace.calcAtomDistance(newangle.i, newangle.j); //angle length
						newangle.theta0 = wspace.calcAtomAngle(newangle.i, newangle.a, newangle.j); //angle angle

						//printf("D %d %lf \n",angle.size(),newangle.k * PhysicsConst::J2kcal * PhysicsConst::Na);
					}
				}else{
						newangle.k = 0.03 * PhysicsConst::kcal2J / PhysicsConst::Na; //force constant per deg
						newangle.l = wspace.calcAtomDistance(newangle.i, newangle.j); //angle length
						newangle.theta0 = wspace.calcAtomAngle(newangle.i, newangle.a, newangle.j); //angle angle
				}
				angle.push_back(newangle);
			}
		}

		if(totalangles != angle.size()) {
			printf("CODE ERROR: TotAng %d != Ang %d in assembleAngleList(); check code\n", totalangles, angle.size());
			return -1;
		}

		for( size_t i = 0; i < wspace.atom.size(); i++) 
		{ 
			// no re collect the equilibrum distances into cov13atom
			for( size_t c = 0; c < atomparam[i].cov13atom.size(); c++) 
			{
				for( aindex = 0; aindex < angle.size(); aindex++) 
				{ 
					// find that bond in anglelist
					if(((angle[aindex].i == i) &&
						(angle[aindex].a == atomparam[i].cov13atom[c].i2) &&
						(angle[aindex].j == atomparam[i].cov13atom[c].i)) ||
						((angle[aindex].j == i) &&
						(angle[aindex].a == atomparam[i].cov13atom[c].i2) && (angle[aindex].i == atomparam[i].cov13atom[c].i)))
						break; // finish for loop !
				}
				if(aindex >= angle.size()) {
					printf("CODE ERROR #2: in assembleAngleList() ; check code\n");
					return -1;
				}


				double db1, db2;
				int ci;

				//find the two bonds involved, firstly a-i
				for(ci = 0; ci < atomparam[angle[aindex].a].cov12atom.size(); ci++) {
					if(atomparam[angle[aindex].a].cov12atom[ci].i == angle[aindex].i)
						break;
				}
				if(ci >= atomparam[angle[aindex].a].cov12atom.size()) {
					printf("CODE ERROR #3a: in assembleAngleList() ; check code\n");
					return -1;
				}
				db1 = (double) atomparam[angle[aindex].a].cov12atom[ci].d;


				//now a-j
				for(ci = 0; ci < atomparam[angle[aindex].a].cov12atom.size(); ci++) {
					if(atomparam[angle[aindex].a].cov12atom[ci].i == angle[aindex].j)
						break;
				}
				if(ci >= atomparam[angle[aindex].a].cov12atom.size()) {
					printf("CODE ERROR #3b: in assembleAngleList() ; check code\n");
					return -1;
				}
				db2 = (double) atomparam[angle[aindex].a].cov12atom[ci].d;


				// now apply the cosine rule to find the equilibrium i-j distance of the angle
				// 'aindex' and save it in the cov13atom array;
				atomparam[i].cov13atom[c].d = sqrt(sqr(db1) + sqr(db2) - 2 * db1 * db2 * cos(angle[aindex].theta0));
				angle[aindex].l = atomparam[i].cov13atom[c].d;
			}
		}


		return 0;
	}


	// Assembles torsion list - impropers are also included in this list from the already
	// established improper[] array
	int FF_Bonded::assembleTorsionList()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		const ParticleStore& atomparam = wspace.atom;
		SnapShotAtom *atom = wspace.cur.atom;

		int i, tindex, c;
		Torsion newtorsion;

		torsion.clear(); // clear previous array

		if(!wspace.ffps().autotorsions)
			return 0;

		size_t totaltorsions = 0;

		for(size_t i = 0; i < wspace.atom.size(); i++) //count total number of torsions
			totaltorsions += atomparam[i].cov14atom.size();

		if((totaltorsions % 2) != 0) {
			THROW(CodeException,"ERROR: Total number of covalent 1,4 links not even - suspect code error/def. file error!");
		}
		totaltorsions /= 2; //half total as every torsion is only recorded once!

		for(i = 0; i < wspace.atom.size(); i++) {
			for(c = 0; c < atomparam[i].cov14atom.size(); c++) { //collect all the covalent links, assume they're all torsion...
				if(i < atomparam[i].cov14atom[c].i)
					continue; //only record half (i.e. don't double record 1-x-x-4 and 4-x-x-1)

				// Found a candidate torsion
				newtorsion.i = i; //record atoms
				newtorsion.a = atomparam[i].cov14atom[c].i3;
				newtorsion.b = atomparam[i].cov14atom[c].i2;
				newtorsion.j = atomparam[i].cov14atom[c].i;
				newtorsion.Type = 0; // indicate that this is a torsion dihedral

				if(!IgnoreForcefieldParams){

					tindex = wspace.ffps().findTorsionType(
						atomparam[newtorsion.i].FFType,
						atomparam[newtorsion.a].FFType,
						atomparam[newtorsion.b].FFType, 
						atomparam[newtorsion.j].FFType);

					if(tindex < 0) 
					{ 
						// ignore if no appropriate torsion Type known
						printf("WARNING: Unknown Torsion Type required in amino acid %d %s (%s %s %s %s) - force field incomplete\n",
							atomparam[newtorsion.i].ir,
							atomparam[newtorsion.i].parentname.c_str(),
							atomparam[newtorsion.i].pdbname.c_str(),
							atomparam[newtorsion.a].pdbname.c_str(),
							atomparam[newtorsion.b].pdbname.c_str(),
							atomparam[newtorsion.j].pdbname.c_str()
							);

						continue; // ignore this torsion, i.e. dont add to list
					} 
					else 
					{
						// now load torsion
						newtorsion.terms = wspace.ffps().TorsionType[tindex].terms;

						for(int iterm = 0; iterm < 4; iterm++) {
							newtorsion.Vn[iterm] = wspace.ffps().TorsionType[tindex].Vn[iterm];
							newtorsion.n[iterm] = wspace.ffps().TorsionType[tindex].n[iterm];
							newtorsion.gamma[iterm] = wspace.ffps().TorsionType[tindex].gamma[iterm];
						}
					}
				}else{
					// just put in 0 parameters
						newtorsion.terms = 0;
						for(int iterm = 0; iterm < 4; iterm++) {
							newtorsion.Vn[iterm] = 0;
							newtorsion.n[iterm] = 2;
							newtorsion.gamma[iterm] = 0;
						}

				}


				torsion.push_back(newtorsion);
			}

		}
		return 0;
	}


	int FF_Bonded::assembleImproperList()
	{
        // Initialisation
        improper.clear();
		if(IgnoreForcefieldParams)
        { 
            return 0;
        }

        int errorCount = 0;

		// Proxies
		WorkSpace& wspace = getWSpace();
		const ParticleStore& atomparam = wspace.atom;
		SnapShotAtom *atom = wspace.cur.atom;

		for(int ir = 0; ir < wspace.res.size(); ir++)
        {
			// Add impropers, (always explicitly defined in residue definition)   
            const MoleculeDefinition* molParam = wspace.res[ir].param;
            const std::vector<DihedralDefinition>& diheds = molParam->improper;

			for(int iimp = 0; iimp < diheds.size(); iimp++) 
            {
                const DihedralDefinition& dihed = diheds[iimp];
				Torsion newimp;

                // If an atom is not found and no backLink\frwdLink is present, ignore the improper...
                if( dihed.roi < 0 ||dihed.roa < 0 || dihed.rob < 0 || dihed.roj < 0 )
                {
                    if(!molParam->hasBackLink() || !molParam->hasFrwdLink())
                    {
                        // DelConnect was previously called to remove the link
                        continue;
                    }
                }

				// printf("Adding improper dihedral..\n");                
				newimp.i = wspace.findParticleBy_ffname(ir + dihed.roi, dihed.ani);
				newimp.a = wspace.findParticleBy_ffname(ir + dihed.roa, dihed.ana);
				newimp.b = wspace.findParticleBy_ffname(ir + dihed.rob, dihed.anb);
				newimp.j = wspace.findParticleBy_ffname(ir + dihed.roj, dihed.anj);

				// check all particles have been found
				if(newimp.i < 0)
                {
					printf("FORCEFIELD ERROR: Improper definition refers to unknown atom res:atnam %d:%s\n",
						ir + dihed.roi, 
                        &dihed.ani[0]
                        );
                    errorCount++;
					continue;
				}
				else if(newimp.a < 0) 
                {
					printf("FORCEFIELD ERROR: Improper definition refers to unknown atom res:atnam %d:%s\n",
						ir + dihed.roa, 
                        &dihed.ana[0]
                        );
                    errorCount++;
					continue;
				}
				else if(newimp.b < 0) 
                {
					printf("FORCEFIELD ERROR: Improper definition refers to unknown atom res:atnam %d:%s\n",
						ir + dihed.rob, 
                        &dihed.anb[0]
                        );
                    errorCount++;
					continue;
				}
				else if(newimp.j < 0) 
                {
					printf("FORCEFIELD ERROR: Improper definition refers to unknown atom res:atnam %d:%s\n",
						ir + dihed.roj, &dihed.anj[0]
                        );
                    errorCount++;
					continue;
				}

				newimp.Type = 1; //mark it's an improper rather than a true torsion (cosmetic)

				int tindex = wspace.ffps().findImproperType(
					atomparam[newimp.i].FFType,
					atomparam[newimp.a].FFType,
					atomparam[newimp.b].FFType, 
					atomparam[newimp.j].FFType);

				if(tindex < 0) 
                { 
                    // ignore if no appropriate torsion Type known
					printf("WARNING: Unknown Improper Type (res %s)- ignoring, forcefield incomplete \n",
						atomparam[newimp.i].parentl3name.c_str());
                    errorCount++;
					continue;
				} 
                else 
                {
					// now load improper (note that in ffps torsions & impropers are treated together so
					// wspace.ffps().TorsionType[tindex]... appears here !
					newimp.terms = wspace.ffps().TorsionType[tindex].terms;
					for(int iterm = 0; iterm < 4; iterm++)
                    {
						newimp.Vn[iterm] = wspace.ffps().TorsionType[tindex].Vn[iterm];
						newimp.n[iterm] = wspace.ffps().TorsionType[tindex].n[iterm];
						newimp.gamma[iterm] = wspace.ffps().TorsionType[tindex].gamma[iterm];
					}
				}
				// newimp now ready - add it to the dihedral list

				improper.push_back(newimp);
			}
		}

		return errorCount;
	}

	Physics::Bond FF_Bonded::getBond( int i, int j ) const
	{
		int index = findBond( i, j );
		if( index == -1 )
		{
			StringBuilder sb;
			sb.setFormat( "The bonded forcefield states that argument atoms %d and %d are not bonded.")(i)(j);
			THROW( ArgumentException, sb.toString() );
		}
		return bond[index];
	}

	bool FF_Bonded::isBonded( int i, int j ) const
	{
		for( size_t k = 0; k < bond.size(); k++ )
		{
			if( ( bond[k].i == i && bond[k].j == j ) || 
				( bond[k].j == i && bond[k].i == j ) )
			{					
				return true;
			}
		}			
		return false;
	}

	void FF_Bonded::resetEquilibriumBonds()
	{
		ASSERT( &getWSpace() != NULL, CodeException, "NULL wspace found in Forcefield::calcForces()");
		SnapShotAtom *atom = getWSpace().cur.atom;

		int                i, j;        // i & j are the atom indices
		int                ib;          // ib is the bond counting variable
		double             d;           // potential energy of the individual bond
		dvector vb;                     // bond dvector  j-->i

		// calculate bonds
		for(ib = 0; ib < bond.size(); ib++) 
		{
			i = bond[ib].i;      // get the atom indices of the bond atoms
			j = bond[ib].j;
			vb.setTo(atom[i].p); // create the bond dvector  j-->i
			vb.sub(atom[j].p);
			d = vb.mag();        // get distance
			bond[ib].l = d; 
		}

	}

	void FF_Bonded::resetEquilibriumAngles()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		SnapShotAtom *atom = wspace.cur.atom;

		int a, i, j, ib;
		dvector iv, jv, uiv, ujv; // bond vectors
		double iv_mag, jv_mag; // and their magnitudes
		double inv_iv_mag, inv_jv_mag;
		double cos_theta;


		for(ib = 0; ib < angle.size(); ib++) {
			i = angle[ib].i;
			a = angle[ib].a;
			j = angle[ib].j;

			iv.x = atom[i].p.x - atom[a].p.x;
			iv.y = atom[i].p.y - atom[a].p.y;
			iv.z = atom[i].p.z - atom[a].p.z;
			iv_mag = sqrt(sqr(iv.x) + sqr(iv.y) + sqr(iv.z));
			inv_iv_mag = 1 / iv_mag;

			uiv.x = iv.x * inv_iv_mag;
			uiv.y = iv.y * inv_iv_mag;
			uiv.z = iv.z * inv_iv_mag;

			jv.x = atom[j].p.x - atom[a].p.x;
			jv.y = atom[j].p.y - atom[a].p.y;
			jv.z = atom[j].p.z - atom[a].p.z;

			jv_mag = sqrt(sqr(jv.x) + sqr(jv.y) + sqr(jv.z));
			inv_jv_mag = 1 / jv_mag;

			ujv.x = jv.x * inv_jv_mag;
			ujv.y = jv.y * inv_jv_mag;
			ujv.z = jv.z * inv_jv_mag;

			cos_theta = uiv.x * ujv.x + uiv.y * ujv.y + uiv.z * ujv.z;

			angle[ib].theta = acos(cos_theta);
			angle[ib].theta0 = angle[ib].theta;
		}
	}

	void FF_Bonded::resetEquilibriumHarmonicDihedrals(){
			calcTorsionForces();
			calcImproperForces();
		

			for( int i = 0; i < torsion.size(); i++) {
				for(int j = 0; j < torsion[i].terms; j++) {
		 			if(torsion[i].n[j]  == 0 ) // set harmonic terms to current phi
					{
						torsion[i].gamma[j] = torsion[i].phi;
					}
				}
			}

			for(int i = 0; i < improper.size() ; i++) { // for every torsion
				for(int j = 0; j < improper[i].terms; j++) {
		 			if(improper[i].n[j]  == 0 ) // set harmonic terms to current phi
					{
						improper[i].gamma[j] = improper[i].phi;
					}
				}	
			}

	}

	int FF_Bonded::findBond(int i, int j) const
	{
		for(int g = 0; g < bond.size(); g++) 
		{
			if((i == bond[g].i) && (j == bond[g].j))
			{
				return g;
			}
			if((i == bond[g].j) && (j == bond[g].i))
			{
				return g;
			}
		}
		return -1;
	}


	int FF_Bonded::findAngle(int i, int a, int j) const
	{
		for(int g = 0; g < angle.size(); g++) 
		{
			if((i == angle[g].i) && (a == angle[g].a) && (j == angle[g].j))
				return g;
			if((i == angle[g].j) && (a == angle[g].a) && (j == angle[g].i))
				return g;
		}
		return -1;
	}


	int FF_Bonded::findTorsion(int i, int a, int b, int j) const
	{
		for(int g = 0; g < torsion.size(); g++) 
		{
			if((i == torsion[g].i) &&
				(a == torsion[g].a) &&
				(b == torsion[g].b) &&
				(j == torsion[g].j)) return g;
			if((i == torsion[g].j) &&
				(a == torsion[g].b) &&
				(b == torsion[g].a) &&
				(j == torsion[g].i)) return g;
		}
		return -1;
	}

	int FF_Bonded::findImproper(int i, int a, int b, int j) const
	{
		for(int g = 0; g < improper.size(); g++) 
		{
			if((i == improper[g].i) &&
				(a == improper[g].a) &&
				(b == improper[g].b) &&
				(j == improper[g].j))
				return g;
		}
		return -1;
	}


	// print Bonded Parameter

	void FF_Bonded::printBondedParameters()
	{
		printBondedParameters(stdout);
	}

	int FF_Bonded::printBondedParameters(char *targetfile)
	{
		FILE *target;
		target = fopen(targetfile, "w");
		printBondedParameters(target);
		fclose(target);
		return 0;
	}

	// -----------------------------------------------
	// SimpleAtomicPrint()
	// DEPRECATED - LEGACY
	// USE class PDB_Writer
	// -----------------------------------------------

	// Kept for compatibility and simplicity.
	int SimpleAtomicPrint(FILE * target, WorkSpace* _wspace)
	{
		char nbuf[20];
		char molchain;
		int alen;
		double zero = 0.0;

		const ParticleStore &atom = _wspace->atom;
		size_t natom = _wspace->atom.size();
		size_t nmolecules = _wspace->mol.size();

		for(size_t i = 0; i < natom; i++) {

			nbuf[0] = ' ';
			if((strlen(atom[i].pdbname.c_str()) >= 4) ||
				((atom[i].pdbname.c_str()[0] >= '1') &&
				 (atom[i].pdbname.c_str()[0] <= '9')))
				strcpy(&nbuf[0], atom[i].pdbname.c_str());
			else
				strcpy(&nbuf[1], atom[i].pdbname.c_str());

			alen = (int) strlen(atom[i].parentl3name.c_str()) - 3;
			if(alen < 0)
				alen = 0;

			if(nmolecules <= 1)
				molchain = ' ';
			else
				molchain = (int)'A' + atom[i].imol;

			fprintf(target, "ATOM  %5d %-4s%c%3s %c%4i%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf",
				i + 1,
				&nbuf[0], ' ',
				&atom[i].parentl3name.c_str()[alen], molchain,
				atom[i].ir, ' ', atom[i].pos().x, atom[i].pos().y, atom[i].pos().z, zero, zero);
			fprintf(target, "\n");
		}

		fprintf(target, "TER\n");
		return 0;
	}

	void FF_Bonded::printBondedParameters(FILE * target)
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		const ParticleStore& atomparam = wspace.atom;
		SnapShotAtom *atom = wspace.cur.atom;
		int i, j;

		fprintf(target, "----------------------------------------------------\n");
		fprintf(target, "Total covalent structure components inc. parameters\n");

		SimpleAtomicPrint(target,&wspace);

		for(i = 0; i < bond.size(); i++) {
			fprintf(target, "Bond %3d: %s <%3d(%-3s) %3d(%-3s)> \tk: %4.2lf (kcal/mol)\t%4.2lf nrml1\n",
				i,
				atomparam[i].parentl3name.c_str(), // ResNam of i
				bond[i].i, atomparam[bond[i].i].pdbname.c_str(), // i, and atname of i
				bond[i].j, atomparam[bond[i].j].pdbname.c_str(), // j, and atname of j
				(double)bond[i].k * PhysicsConst::J2kcal * PhysicsConst::Na,(double) bond[i].l); // and bond length

		}
		fprintf(target, " Total bonds: %d\n", bond.size());

		for(i = 0; i < angle.size(); i++) {
			fprintf(target, "Angle %3d: %s <%3d(%-3s) %3d(%-3s) %3d(%-3s)> k: %5.2lf (kcal/mol/rad) l:%4.1lfA %5.1lfdeg\n",
				i,
				atomparam[i].parentl3name.c_str(), // ResNam of i
				angle[i].i, atomparam[angle[i].i].pdbname.c_str(), // i, and atname of i
				angle[i].a, atomparam[angle[i].a].pdbname.c_str(), // a, and atname of a
				angle[i].j, atomparam[angle[i].j].pdbname.c_str(), // j, and atname of j
				angle[i].k * PhysicsConst::J2kcal * PhysicsConst::Na, angle[i].l, angle[i].theta0 * 180 / Maths::MathConst::PI);
		}
		fprintf(target, " Total angles: %d\n", angle.size());

		for(i = 0; i < torsion.size(); i++) {
			for(j = 0; j < torsion[i].terms; j++) {
				if(j == 0)
					fprintf(target, "Torsion %3d: %3s <%3d(%-3s) %3d(%-3s) %3d(%-3s) %3d(%-3s)> ",
					i,
					atomparam[i].parentl3name.c_str(), // ResNam of i
					torsion[i].i, atomparam[torsion[i].i].pdbname.c_str(), // i, and atname of i
					torsion[i].a, atomparam[torsion[i].a].pdbname.c_str(), // a, and atname of a
					torsion[i].a, atomparam[torsion[i].b].pdbname.c_str(), // b, and atname of b
					torsion[i].j, atomparam[torsion[i].j].pdbname.c_str()); // j, and atname of j
				else
					fprintf(target, " ");
				fprintf(target, " Vn: %5.2lf n:%2.lf gamma: %.1lf \n",
					torsion[i].Vn[j] * PhysicsConst::J2kcal * PhysicsConst::Na, torsion[i].n[j], torsion[i].gamma[j] * 180.0 / Maths::MathConst::PI);

			}
		}
		fprintf(target, " Total torsions: %d\n", torsion.size());

		for(i=0;i<improper.size();i++) {
			for(j=0;j<improper[i].terms;j++){
				if(j == 0)
					fprintf(target, "Improper %3d: %3s <%3d(%-3s) %3d(%-3s) %3d(%-3s) %3d(%-3s)> ",
					i,
					atomparam[i].parentl3name.c_str(), // ResNam of i
					improper[i].i, atomparam[improper[i].i].pdbname.c_str(), // i, and atname of i
					improper[i].a, atomparam[improper[i].a].pdbname.c_str(), // a, and atname of a
					improper[i].a, atomparam[improper[i].b].pdbname.c_str(), // b, and atname of b
					improper[i].j, atomparam[improper[i].j].pdbname.c_str()); // j, and atname of j
				else fprintf(target," ");
				fprintf(target," Vn: %5.2lf n:%2.lf gamma: %.1lf \n",
					improper[i].Vn[j]*PhysicsConst::J2kcal*PhysicsConst::Na,
					improper[i].n[j],
					improper[i].gamma[j]*180/Maths::MathConst::PI);
			}
		}
		fprintf(target," Total impropers: %d\n",improper.size());
	}


	void FF_Bonded::calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level)
	{
		if(level > Summary) 
		{
			if(DoBonds)
				calcBondEnergies_Verbose();
			if(DoAngles)
				calcAngleEnergies_Verbose();
			if(DoTorsions)
				calcTorsionEnergies_Verbose();
			if(DoImpropers)
				calcImproperEnergies_Verbose();
		} 
		else 
		{
			WorkSpace& wspace = getWSpace();
			calcEnergies();
			printf(" Bonds:         %10.3lf kcal/mol\n", double(wspace.ene.epot_bond) * PhysicsConst::J2kcal * PhysicsConst::Na);
			printf(" Angles:        %10.3lf kcal/mol\n", double(wspace.ene.epot_angle) * PhysicsConst::J2kcal * PhysicsConst::Na);
			printf(" Tors/Impr:     %10.3lf kcal/mol\n", double(wspace.ene.epot_torsion) * PhysicsConst::J2kcal * PhysicsConst::Na);
		}
	}

	void FF_Bonded::calcEnergies()
	{
		calcForces();
	}

	void FF_Bonded::calcForces()
	{
		if(DoBonds)
			calcBondForces();
		if(DoAngles)
			calcAngleForces();
		if(DoTorsions)
			calcTorsionForces();
		if(DoImpropers)
			calcImproperForces();

		epot = epot_bond +
			epot_angle +
			epot_torsion +
			epot_improper;
	}

	void FF_Bonded::info() const
	{ 
		// prints a little block of parameter information
		ForcefieldBase::info(); // standard info 
		if(DoBonds)
			printf(" Bonds:     %8d\n", bond.size());
		if(DoAngles)
			printf(" Angles:    %8d\n", angle.size());
		if(DoTorsions)
			printf(" Torsions:  %8d\n", torsion.size());
		if(DoImpropers)
			printf(" Impropers: %8d\n", improper.size());

		printf("\n");
	}

	void FF_Bonded::infoLine() const 
	{ 
		// prints a line of current energies
		const WorkSpace& wspace = getWSpace();
		const Hamiltonian *printene = &wspace.ene;
		if(DoBonds)
			printf("\t%6.1lf", printene->epot_bond * PhysicsConst::J2kcal * PhysicsConst::Na);
		if(DoAngles)
			printf("\t%6.1lf", printene->epot_angle * PhysicsConst::J2kcal * PhysicsConst::Na);
		if(DoTorsions)
			printf("\t%6.1lf", printene->epot_torsion * PhysicsConst::J2kcal * PhysicsConst::Na);
	}

	void FF_Bonded::infoLineHeader() const 
	{
		// prints the headers for the above function
		if(DoBonds)
			printf("\t%6s", "Ebond");
		if(DoAngles)
			printf("\t%6s", "Eangle");
		if(DoTorsions)
			printf("\t%6s", "Edihed");
	}

	// -------------------------------------------------------------------------------
	// Class: FF_Bonded
	// Function: calcBondForces();
	// -------------------------------------------------------------------------------
	// Parameters: none
	//
	// from "wspace.h" :
	// typedef struct __Bond{ // bond i-j
	// int i,j; // indices of atoms
	// double l; // bond length
	// double k; // force constant
	// double d; // actual length of bond
	// } Bond;

	// Handles harmonic bond stretching
	// It loops through all bonds defined in the array newbond and calculates
	// the potential energy as well as the force acting along the i-j dvector 
	// The forces & potentials are then added to the atoms

	void FF_Bonded::calcBondForces()
	{
		ASSERT( &getWSpace() != NULL, CodeException, "NULL wspace found in Forcefield::calcForces()");

		// Proxies
		WorkSpace& wspace = getWSpace();
		ParticleStore& atomparam = wspace.atom;
		SnapShotAtom *atom = wspace.cur.atom;

		int                i, j;        // i & j are the atom indices
		int                ib;          // ib is the bond counting variable
		double             force_magnitude, k;       // force magnitude, force constant
		dvector vb;          // bond dvector  j-->i
		double             d,potential; // potential energy of the individual bond

		epot_bond = 0;

		// calculate bonds
		for(ib = 0; ib < bond.size(); ib++) 
		{
			if(Scope==OnlyBackbone) 
			{
				if(!wspace.atom[bond[ib].i].isBackbone())
					continue;
				if(!wspace.atom[bond[ib].j].isBackbone())
					continue;
			}

			i = bond[ib].i; // get the atom indices of the bond atoms
			j = bond[ib].j;
			k = bond[ib].k; // get force constant in J nrml1-2

			vb.setTo(atom[i].p); // create the bond dvector  j-->i
			vb.sub(atom[j].p);
			d = vb.mag();        // get distance

			force_magnitude = 2 * k * (d - bond[ib].l) / PhysicsConst::Angstrom; //calculate force
		
			// Virial contribution
			wspace.ene.InternalVirial += d * force_magnitude* PhysicsConst::Angstrom;

			potential = k * sqr(d - bond[ib].l);
			epot_bond += potential;

			vb.mul(force_magnitude / d); // divide by length (unify) and times by force_magnitude(forcemagnitude)
			

			atom[i].f.sub(vb); // add force and reaction force to the two bond atoms
			atom[j].f.add(vb);
			atomparam[i].epot += potential * 0.5; // add half of potential to each atom
			atomparam[j].epot += potential * 0.5;
		}

		wspace.ene.epot_bond += epot_bond; // add up total potential of system
		wspace.ene.epot += epot_bond;
	}


	void FF_Bonded::calcBondEnergies_Verbose()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		const ParticleStore& atomparam = wspace.atom;
		SnapShotAtom *atom = wspace.cur.atom;

		int    i, j, ib;
		double k,d; //force magnitude, force constant
		double potential;
		double potential_sum = 0;

		wspace.ene.epot_bond = 0;

		// calculate bonds
		printf("\n i Nr(Name) j Nr(Name) Length Eq.Lenh k V (kcal/mol)\n");
		for(ib = 0; ib < bond.size(); ib++) {
			i = bond[ib].i;
			j = bond[ib].j;
			d = dist(atom[i].p,atom[j].p);
			k = bond[ib].k; //force constant in J nrml1-2

			potential = k * sqr(d - bond[ib].l);
			potential_sum += potential;
			wspace.ene.epot_bond += potential;

			printf("Bond: %5d(%4s)%5d(%4s) %5.2lf %5.2lf %8.3lf %8.3lf\n",
				i, atomparam[i].pdbname.c_str(),
				j, atomparam[j].pdbname.c_str(),
				d,
				bond[ib].l,
				k * PhysicsConst::J2kcal * PhysicsConst::Na,
				potential * PhysicsConst::J2kcal * PhysicsConst::Na);

		}
		printf(" ----------------\n");
		printf(" Total: %8.3lf \n\n", potential_sum * PhysicsConst::J2kcal * PhysicsConst::Na);
		wspace.ene.epot += wspace.ene.epot_bond;
	}



	// -------------------------------------------------------------------------------
	// Class: FF_Bonded
	// Function: calcAngleForces();
	// -------------------------------------------------------------------------------
	// Parameters: none
	//
	// As above but streamlined, the equivlent slow instructions are commented out
	//

	void FF_Bonded::calcAngleForces()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		ParticleStore& atomparam = wspace.atom;
		SnapShotAtom *atom = wspace.cur.atom;

		int a, i, j, ib;
		double k, potential, force_magnitude;
		dvector iv, jv, uiv, ujv; // bond vectors
		double iv_mag, jv_mag; // and their magnitudes
		double inv_iv_mag, inv_jv_mag;
		double sin_theta, cos_theta;
		dvector forcev; // force dvector 
		dvector bondnormal; // normal dvector  (plane of bond)

		epot_angle = 0;

		for(ib = 0; ib < angle.size(); ib++) 
		{
			if( angle[ib].k == 0 ) continue;
			if(Scope==OnlyBackbone) 
			{
				if(!wspace.atom[angle[ib].i].isBackbone())
					continue;
				if(!wspace.atom[angle[ib].a].isBackbone())
					continue;
				if(!wspace.atom[angle[ib].j].isBackbone())
					continue;
			}

			i = angle[ib].i;
			a = angle[ib].a;
			j = angle[ib].j;

			//iv.setTo(atom[i].p);
			//iv.sub(atom[a].p);
			iv.x = atom[i].p.x - atom[a].p.x;
			iv.y = atom[i].p.y - atom[a].p.y;
			iv.z = atom[i].p.z - atom[a].p.z;
			//iv_mag = iv.mag();
			iv_mag = sqrt(sqr(iv.x) + sqr(iv.y) + sqr(iv.z));
			inv_iv_mag = 1 / iv_mag;

			//uiv.setTo(iv);
			//uiv.div(iv_mag);
			uiv.x = iv.x * inv_iv_mag;
			uiv.y = iv.y * inv_iv_mag;
			uiv.z = iv.z * inv_iv_mag;

			//jv.setTo(atom[j].p);
			//jv.sub(atom[a].p);
			jv.x = atom[j].p.x - atom[a].p.x;
			jv.y = atom[j].p.y - atom[a].p.y;
			jv.z = atom[j].p.z - atom[a].p.z;

			//jv_mag = jv.mag();
			jv_mag = sqrt(sqr(jv.x) + sqr(jv.y) + sqr(jv.z));
			inv_jv_mag = 1 / jv_mag;

			//ujv.setTo(jv);
			//ujv.div(jv_mag);
			ujv.x = jv.x * inv_jv_mag;
			ujv.y = jv.y * inv_jv_mag;
			ujv.z = jv.z * inv_jv_mag;

			//cos_theta = uiv.scalarProduct(ujv);
			cos_theta = uiv.x * ujv.x + uiv.y * ujv.y + uiv.z * ujv.z;

			angle[ib].theta = acos(cos_theta);
			sin_theta = sin(angle[ib].theta);

			k = angle[ib].k;
			potential = k * sqr(angle[ib].theta - angle[ib].theta0);
			epot_angle += potential;

			force_magnitude = 2 * k * (angle[ib].theta - angle[ib].theta0) / (PhysicsConst::Angstrom);

			//bondnormal.setTo(uiv);
			//bondnormal.mul(cos_theta);
			//bondnormal.sub(ujv);
			//bondnormal.div(sin_theta);
			//bondnormal.mul(force_magnitude/iv_mag);

			bondnormal.x = force_magnitude * inv_iv_mag * (uiv.x * cos_theta - ujv.x) / sin_theta;
			bondnormal.y = force_magnitude * inv_iv_mag * (uiv.y * cos_theta - ujv.y) / sin_theta;
			bondnormal.z = force_magnitude * inv_iv_mag * (uiv.z * cos_theta - ujv.z) / sin_theta;


			//atom[i].f.sub(bondnormal);
			atom[i].f.x -= bondnormal.x;
			atom[i].f.y -= bondnormal.y;
			atom[i].f.z -= bondnormal.z;

			//atom[a].f.add(bondnormal);
			atom[a].f.x += bondnormal.x;
			atom[a].f.y += bondnormal.y;
			atom[a].f.z += bondnormal.z;


			//bondnormal.setTo(ujv);
			//bondnormal.mul(cos_theta);
			//bondnormal.sub(uiv);
			//bondnormal.div(sin_theta);
			//bondnormal.mul(force_magnitude/jv_mag);

			bondnormal.x = force_magnitude * inv_jv_mag * (ujv.x * cos_theta - uiv.x) / sin_theta;
			bondnormal.y = force_magnitude * inv_jv_mag * (ujv.y * cos_theta - uiv.y) / sin_theta;
			bondnormal.z = force_magnitude * inv_jv_mag * (ujv.z * cos_theta - uiv.z) / sin_theta;

			//atom[j].f.sub(bondnormal);
			atom[j].f.x -= bondnormal.x;
			atom[j].f.y -= bondnormal.y;
			atom[j].f.z -= bondnormal.z;

			//atom[a].f.add(bondnormal);
			atom[a].f.x += bondnormal.x;
			atom[a].f.y += bondnormal.y;
			atom[a].f.z += bondnormal.z;

			atomparam[i].epot += potential*0.25;
			atomparam[j].epot += potential*0.25;
			atomparam[a].epot += potential*0.50;
		}

		wspace.ene.epot_angle += epot_angle;
		wspace.ene.epot += epot_angle;
	}








	// calculate angles
	void FF_Bonded::calcAngleEnergies_Verbose()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		const ParticleStore& atomparam = wspace.atom;
		SnapShotAtom *atom = wspace.cur.atom;

		int a, i, j, ib;
		double k, potential;
		double potential_sum = 0;
		dvector iv, jv, uiv, ujv; // bond vectors
		double iv_mag, jv_mag; // and their magnitudes
		double sin_theta, cos_theta;
		dvector forcev; // force dvector 
		dvector bondnormal; // normal dvector  (plane of bond)

		wspace.ene.epot_angle = 0;
		printf("\n i Nr(Name) a nr(Atom) j Nr(Name) theta theta0 k ene (kcal/mol)\n");

		for(ib = 0; ib < angle.size(); ib++) {
			if( angle[ib].k == 0 ) continue;
			i = angle[ib].i;
			a = angle[ib].a;
			j = angle[ib].j;

			iv.setTo(atom[i].p);
			iv.sub(atom[a].p);
			iv_mag = iv.mag();
			uiv.setTo(iv);
			uiv.div(iv_mag);

			jv.setTo(atom[j].p);
			jv.sub(atom[a].p);
			jv_mag = jv.mag();
			ujv.setTo(jv);
			ujv.div(jv_mag);

			cos_theta = uiv.scalarProduct(ujv);
			angle[ib].theta = acos(cos_theta);
			sin_theta = sin(angle[ib].theta);

			k = angle[ib].k;

			potential = k * sqr(angle[ib].theta - angle[ib].theta0);

			wspace.ene.epot_angle += potential;
			potential_sum += potential;

			printf("Angle: %5d(%4s)%5d(%4s)%5d(%4s) %5.2lf %5.2lf %8.3lf %8.3lf\n",
				i, atomparam[i].pdbname.c_str(),
				a, atomparam[a].pdbname.c_str(),
				j, atomparam[j].pdbname.c_str(),
				angle[ib].theta * 180.0 / Maths::MathConst::PI, angle[ib].theta0 * 180.0 / Maths::MathConst::PI,
				k * PhysicsConst::J2kcal * PhysicsConst::Na / sqr(180.0) * sqr(Maths::MathConst::PI), potential * PhysicsConst::J2kcal * PhysicsConst::Na);
		}
		printf(" ----------------\n");
		printf(" Total: %8.3lf\n\n",
			potential_sum * PhysicsConst::J2kcal * PhysicsConst::Na);

		wspace.ene.epot += wspace.ene.epot_angle;
	}



















	/// Torsion forces

	template <bool verbose, bool passive>
	void calcDihedralForces(
		WorkSpace& wspace, 
		Torsion &dihedral, 
		double &epot_dihedral
	)
	{
		// Proxies
		ParticleStore& atomparam = wspace.atom;
		SnapShotAtom *atom = wspace.cur.atom;

		double  epot_dihedral_this=0;
		int     ti, ta, tb, tj;
		double  phi, sin_phi, cos_phi;
		dvector vti_vta, vta_vtb, vtb_vtj;
		dvector nrml1, nrml2, nrml3;
		double  inv_nrml1_mag, inv_nrml2_mag, inv_nrml3_mag;

		dvector dcosdnrml1;
		dvector dcosdnrml2;
		dvector dsindnrml3;
		dvector dsindnrml2;
		dvector f, fi, fab, fj;

		ti = dihedral.i;
		ta = dihedral.a;
		tb = dihedral.b;
		tj = dihedral.j;

		vti_vta.x = atom[ti].p.x - atom[ta].p.x;
		vti_vta.y = atom[ti].p.y - atom[ta].p.y;
		vti_vta.z = atom[ti].p.z - atom[ta].p.z;

		vta_vtb.x = atom[ta].p.x - atom[tb].p.x;
		vta_vtb.y = atom[ta].p.y - atom[tb].p.y;
		vta_vtb.z = atom[ta].p.z - atom[tb].p.z;

		vtb_vtj.x = atom[tb].p.x - atom[tj].p.x;
		vtb_vtj.y = atom[tb].p.y - atom[tj].p.y;
		vtb_vtj.z = atom[tb].p.z - atom[tj].p.z;

		fi.x = 0;
		fi.y = 0;
		fi.z = 0;
		fab.x = 0;
		fab.y = 0;
		fab.z = 0;
		fj.x = 0;
		fj.y = 0;
		fj.z = 0;

		nrml1.x = (vti_vta.y * vta_vtb.z - vti_vta.z * vta_vtb.y);
		nrml1.y = (vti_vta.z * vta_vtb.x - vti_vta.x * vta_vtb.z);
		nrml1.z = (vti_vta.x * vta_vtb.y - vti_vta.y * vta_vtb.x);

		nrml2.x = (vta_vtb.y * vtb_vtj.z - vta_vtb.z * vtb_vtj.y);
		nrml2.y = (vta_vtb.z * vtb_vtj.x - vta_vtb.x * vtb_vtj.z);
		nrml2.z = (vta_vtb.x * vtb_vtj.y - vta_vtb.y * vtb_vtj.x);

		nrml3.x = (vta_vtb.y * nrml1.z - vta_vtb.z * nrml1.y);
		nrml3.y = (vta_vtb.z * nrml1.x - vta_vtb.x * nrml1.z);
		nrml3.z = (vta_vtb.x * nrml1.y - vta_vtb.y * nrml1.x);

		inv_nrml1_mag = 1.0 / sqrt(sqr(nrml1.x) + sqr(nrml1.y) + sqr(nrml1.z));
		inv_nrml2_mag = 1.0 / sqrt(sqr(nrml2.x) + sqr(nrml2.y) + sqr(nrml2.z));
		inv_nrml3_mag = 1.0 / sqrt(sqr(nrml3.x) + sqr(nrml3.y) + sqr(nrml3.z));


		cos_phi = (nrml1.x * nrml2.x + nrml1.y * nrml2.y + nrml1.z * nrml2.z) * inv_nrml1_mag * inv_nrml2_mag;
		sin_phi = (nrml3.x * nrml2.x + nrml3.y * nrml2.y + nrml3.z * nrml2.z) * inv_nrml3_mag * inv_nrml2_mag;

		nrml2.x *= inv_nrml2_mag;
		nrml2.y *= inv_nrml2_mag;
		nrml2.z *= inv_nrml2_mag;

		phi = -atan2(sin_phi, cos_phi);

		if(fabs(sin_phi) > 0.1) {
			nrml1.x *= inv_nrml1_mag;
			nrml1.y *= inv_nrml1_mag;
			nrml1.z *= inv_nrml1_mag;

			dcosdnrml1.x = inv_nrml1_mag * (nrml1.x * cos_phi - nrml2.x);
			dcosdnrml1.y = inv_nrml1_mag * (nrml1.y * cos_phi - nrml2.y);
			dcosdnrml1.z = inv_nrml1_mag * (nrml1.z * cos_phi - nrml2.z);

			dcosdnrml2.x = inv_nrml2_mag * (nrml2.x * cos_phi - nrml1.x);
			dcosdnrml2.y = inv_nrml2_mag * (nrml2.y * cos_phi - nrml1.y);
			dcosdnrml2.z = inv_nrml2_mag * (nrml2.z * cos_phi - nrml1.z);

		} else {
			nrml3.x *= inv_nrml3_mag;
			nrml3.y *= inv_nrml3_mag;
			nrml3.z *= inv_nrml3_mag;

			dsindnrml3.x = inv_nrml3_mag * (nrml3.x * sin_phi - nrml2.x);
			dsindnrml3.y = inv_nrml3_mag * (nrml3.y * sin_phi - nrml2.y);
			dsindnrml3.z = inv_nrml3_mag * (nrml3.z * sin_phi - nrml2.z);

			dsindnrml2.x = inv_nrml2_mag * (nrml2.x * sin_phi - nrml3.x);
			dsindnrml2.y = inv_nrml2_mag * (nrml2.y * sin_phi - nrml3.y);
			dsindnrml2.z = inv_nrml2_mag * (nrml2.z * sin_phi - nrml3.z);
		}

		dihedral.phi = phi;

		for(int j = 0; j < dihedral.terms; j++) {
			double k = dihedral.Vn[j];
			double n = dihedral.n[j];
			double gamma = dihedral.gamma[j];
			double epot, dEdphi;

			if(n > 0) { // sin potential
				epot = k * (1.0 + cos(n * phi + gamma));
				dEdphi = -n * k * sin(n * phi + gamma);
			} else { // harmonic potential
				double diff = phi - gamma;
				if(diff < -Maths::MathConst::PI)
					diff += 2.0 * Maths::MathConst::PI;
				else if(diff > Maths::MathConst::PI)
					diff -= 2.0 * Maths::MathConst::PI;
				epot = k * sqr(diff);
				dEdphi = 2.0 * k * diff;
			}
			epot_dihedral_this += epot;
			atomparam[ti].epot += epot*0.25;
			atomparam[ta].epot += epot*0.25;
			atomparam[tb].epot += epot*0.25;
			atomparam[tj].epot += epot*0.25;

			// forces
			if(fabs(sin_phi) > 0.1) {
				dEdphi /= sin_phi;
				fi.x += dEdphi * (vta_vtb.y * dcosdnrml1.z - vta_vtb.z * dcosdnrml1.y);
				fi.y += dEdphi * (vta_vtb.z * dcosdnrml1.x - vta_vtb.x * dcosdnrml1.z);
				fi.z += dEdphi * (vta_vtb.x * dcosdnrml1.y - vta_vtb.y * dcosdnrml1.x);

				fj.x += dEdphi * (vta_vtb.z * dcosdnrml2.y - vta_vtb.y * dcosdnrml2.z);
				fj.y += dEdphi * (vta_vtb.x * dcosdnrml2.z - vta_vtb.z * dcosdnrml2.x);
				fj.z += dEdphi * (vta_vtb.y * dcosdnrml2.x - vta_vtb.x * dcosdnrml2.y);

				fab.x += dEdphi * (vti_vta.z * dcosdnrml1.y - vti_vta.y * dcosdnrml1.z + vtb_vtj.y * dcosdnrml2.z - vtb_vtj.z * dcosdnrml2.y);
				fab.y += dEdphi * (vti_vta.x * dcosdnrml1.z - vti_vta.z * dcosdnrml1.x + vtb_vtj.z * dcosdnrml2.x - vtb_vtj.x * dcosdnrml2.z);
				fab.z += dEdphi * (vti_vta.y * dcosdnrml1.x - vti_vta.x * dcosdnrml1.y + vtb_vtj.x * dcosdnrml2.y - vtb_vtj.y * dcosdnrml2.x);

			} else {
				dEdphi /= -cos_phi;

				fi.x += dEdphi * ((vta_vtb.y * vta_vtb.y + vta_vtb.z * vta_vtb.z) * dsindnrml3.x - vta_vtb.x * vta_vtb.y * dsindnrml3.y - vta_vtb.x * vta_vtb.z * dsindnrml3.z);
				fi.y += dEdphi * ((vta_vtb.z * vta_vtb.z + vta_vtb.x * vta_vtb.x) * dsindnrml3.y - vta_vtb.y * vta_vtb.z * dsindnrml3.z - vta_vtb.y * vta_vtb.x * dsindnrml3.x);
				fi.z += dEdphi * ((vta_vtb.x * vta_vtb.x + vta_vtb.y * vta_vtb.y) * dsindnrml3.z - vta_vtb.z * vta_vtb.x * dsindnrml3.x - vta_vtb.z * vta_vtb.y * dsindnrml3.y);

				fj.x += dEdphi * (dsindnrml2.y * vta_vtb.z - dsindnrml2.z * vta_vtb.y);
				fj.y += dEdphi * (dsindnrml2.z * vta_vtb.x - dsindnrml2.x * vta_vtb.z);
				fj.z += dEdphi * (dsindnrml2.x * vta_vtb.y - dsindnrml2.y * vta_vtb.x);

				fab.x += dEdphi * (-(vta_vtb.y * vti_vta.y + vta_vtb.z * vti_vta.z) * dsindnrml3.x
					+ (2.0 * vta_vtb.x * vti_vta.y - vti_vta.x * vta_vtb.y) * dsindnrml3.y
					+ (2.0 * vta_vtb.x * vti_vta.z - vti_vta.x * vta_vtb.z) * dsindnrml3.z + dsindnrml2.z * vtb_vtj.y - dsindnrml2.y * vtb_vtj.z);
				fab.y += dEdphi * (-(vta_vtb.z * vti_vta.z + vta_vtb.x * vti_vta.x) * dsindnrml3.y
					+ (2.0 * vta_vtb.y * vti_vta.z - vti_vta.y * vta_vtb.z) * dsindnrml3.z
					+ (2.0 * vta_vtb.y * vti_vta.x - vti_vta.y * vta_vtb.x) * dsindnrml3.x + dsindnrml2.x * vtb_vtj.z - dsindnrml2.z * vtb_vtj.x);
				fab.z += dEdphi * (-(vta_vtb.x * vti_vta.x + vta_vtb.y * vti_vta.y) * dsindnrml3.z
					+ (2.0 * vta_vtb.z * vti_vta.x - vti_vta.z * vta_vtb.x) * dsindnrml3.x
					+ (2.0 * vta_vtb.z * vti_vta.y - vti_vta.z * vta_vtb.y) * dsindnrml3.y + dsindnrml2.y * vtb_vtj.x - dsindnrml2.x * vtb_vtj.y);
			}
		} 

		// passive forcefields do not deposit forces or energies
		if(!passive){
			fi.mul(PhysicsConst::invAngstrom);
			fab.mul(PhysicsConst::invAngstrom);
			fj.mul(PhysicsConst::invAngstrom);

			atom[ti].f.x += fi.x;
			atom[ti].f.y += fi.y;
			atom[ti].f.z += fi.z;

			atom[ta].f.x += fab.x - fi.x;
			atom[ta].f.y += fab.y - fi.y;
			atom[ta].f.z += fab.z - fi.z;

			atom[tb].f.x += fj.x - fab.x;
			atom[tb].f.y += fj.y - fab.y;
			atom[tb].f.z += fj.z - fab.z;

			atom[tj].f.x -= fj.x;
			atom[tj].f.y -= fj.y;
			atom[tj].f.z -= fj.z;
		}
		// the energy is returned in any case, the caller must decide wether
		// to add it to the total energy of the workspace

		epot_dihedral += epot_dihedral_this;

		if(verbose){
		  // Look if the dihedral is a peptide bond dihedral
			std::string pepbuf("     ");
			if((strcmp(atomparam[dihedral.i].pdbname.c_str(), "CA") == 0) &&
				 (strcmp(atomparam[dihedral.j].pdbname.c_str(), "CA") == 0)) {

					if((strcmp(atomparam[dihedral.a].pdbname.c_str(), "N") == 0) &&
						(strcmp(atomparam[dihedral.b].pdbname.c_str(), "C") == 0)) {
							pepbuf = "(pep)";
					}
					if((strcmp(atomparam[dihedral.b].pdbname.c_str(), "N") == 0) &&
						(strcmp(atomparam[dihedral.a].pdbname.c_str(), "C") == 0)) {
							pepbuf = "(pep)";
					}
			}


			printf("dihedral: %4s %5d(%4s)[%3s] %5d(%4s)[%3s] %5d(%4s)[%3s] ",
				atomparam[dihedral.a].parentl3name.c_str(),
				dihedral.i, atomparam[dihedral.i].pdbname.c_str(),
				wspace.ffps().AtomType[atomparam[dihedral.i].FFType].name.c_str(),
				dihedral.a, atomparam[dihedral.a].pdbname.c_str(),
				wspace.ffps().AtomType[atomparam[dihedral.a].FFType].name.c_str(),
				dihedral.b, atomparam[dihedral.b].pdbname.c_str(),
				wspace.ffps().AtomType[atomparam[dihedral.b].FFType].name.c_str());
			if( ( dihedral.terms == 1) && (dihedral.n[0] == 0 ) ){   // a simple restraint ??
				printf(" %5d(%4s)[%3s] %s %8.3lf (gamma=%8.3lf) (k=%8.3lf) %8.3lf\n",
					dihedral.j, atomparam[dihedral.j].pdbname.c_str(),
					wspace.ffps().AtomType[atomparam[dihedral.j].FFType].name.c_str(),
					pepbuf.c_str(), 
					phi * 180.0 / Maths::MathConst::PI,
					dihedral.gamma[0] * 180.0 / Maths::MathConst::PI,
					dihedral.Vn[0] * PhysicsConst::J2kcal * PhysicsConst::Na,
					epot_dihedral_this * PhysicsConst::J2kcal * PhysicsConst::Na);
			}else{
				printf(" %5d(%4s)[%3s] %s %8.3lf %8.3lf\n",
					dihedral.j, atomparam[dihedral.j].pdbname.c_str(),
					wspace.ffps().AtomType[atomparam[dihedral.j].FFType].name.c_str(),
					pepbuf.c_str(), phi * 180.0 / Maths::MathConst::PI,
					epot_dihedral_this * PhysicsConst::J2kcal * PhysicsConst::Na);
			}
		}

	}

	void calcDihedralForcesNonVerbose(WorkSpace &wspace, Torsion &dihedral, double &epot_dihedral){
		calcDihedralForces<false,false>(wspace,dihedral,epot_dihedral);
	}
	void calcDihedralForcesVerbose(WorkSpace &wspace, Torsion &dihedral, double &epot_dihedral){
		calcDihedralForces<true,false>(wspace,dihedral,epot_dihedral);
	}
	void calcDihedralForcesNonVerbosePassive(WorkSpace &wspace, Torsion &dihedral, double &epot_dihedral){
		calcDihedralForces<false,true>(wspace,dihedral,epot_dihedral);
	}
	void calcDihedralForcesVerbosePassive(WorkSpace &wspace, Torsion &dihedral, double &epot_dihedral){
		calcDihedralForces<true,true>(wspace,dihedral,epot_dihedral);
	}






	void FF_Bonded::calcTorsionForces()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		epot_torsion = 0;

		for( int i = 0; i < torsion.size(); i++) 
		{ 
			// for every torsion
			if(Scope==OnlyBackbone) 
			{
				if(!wspace.atom[torsion[i].i].isBackbone())
					continue;
				if(!wspace.atom[torsion[i].a].isBackbone())
					continue;
				if(!wspace.atom[torsion[i].b].isBackbone())
					continue;
				if(!wspace.atom[torsion[i].j].isBackbone())
					continue;
			}

			calcDihedralForces<false,false>(wspace,torsion[i],epot_torsion);
		}
		wspace.ene.epot_torsion += epot_torsion;
		wspace.ene.epot += epot_torsion;

	}


	void FF_Bonded::calcTorsionEnergies_Verbose()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		epot_torsion = 0;

		printf("\n i Nr(Name) a nr(Atom) b nr(Atom) j Nr(Name) ene (kcal/mol)\n");
		for( int i = 0; i < torsion.size(); i++) 
		{ 
			// for every torsion
			if(Scope==OnlyBackbone) 
			{
				if(!wspace.atom[torsion[i].i].isBackbone())
					continue;
				if(!wspace.atom[torsion[i].a].isBackbone())
					continue;
				if(!wspace.atom[torsion[i].b].isBackbone())
					continue;
				if(!wspace.atom[torsion[i].j].isBackbone())
					continue;
			}

			calcDihedralForces<true,false>(wspace,torsion[i],epot_torsion);
		}
		printf(" ----------------\n");
		printf(" Total: %8.3lf\n\n", epot_torsion * PhysicsConst::J2kcal * PhysicsConst::Na);

		wspace.ene.epot_torsion += epot_torsion;
		wspace.ene.epot += epot_torsion;
	}


	void FF_Bonded::calcImproperForces()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		epot_improper = 0;

		for(int i = 0; i < improper.size() ; i++) { // for every torsion
			if(Scope==OnlyBackbone) 
			{
				if(!wspace.atom[torsion[i].i].isBackbone())
					continue;
				if(!wspace.atom[torsion[i].a].isBackbone())
					continue;
				if(!wspace.atom[torsion[i].b].isBackbone())
					continue;
				if(!wspace.atom[torsion[i].j].isBackbone())
					continue;
			}

			calcDihedralForces<false,false>(wspace,improper[i],epot_improper);
		}
		wspace.ene.epot_torsion += epot_improper;
		wspace.ene.epot += epot_improper;
	}


	void FF_Bonded::calcImproperEnergies_Verbose()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		epot_improper = 0;

		printf("\n i Nr(Name) a nr(Atom) b nr(Atom) j Nr(Name) ene (kcal/mol)\n");
		for(int i = 0; i < improper.size() ; i++) { // for every torsion
			if(Scope==OnlyBackbone) 
			{
				if(!wspace.atom[torsion[i].i].isBackbone())
					continue;
				if(!wspace.atom[torsion[i].a].isBackbone())
					continue;
				if(!wspace.atom[torsion[i].b].isBackbone())
					continue;
				if(!wspace.atom[torsion[i].j].isBackbone())
					continue;
			}

			calcDihedralForces<true,false>(wspace,improper[i],epot_improper);
		}
		printf(" ----------------\n");
		printf(" Total: %8.3lf\n\n", epot_improper * PhysicsConst::J2kcal * PhysicsConst::Na);

		wspace.ene.epot_torsion += epot_improper;
		wspace.ene.epot += epot_improper;
	}




	int FF_Bonded::printPSFfile_bondedparams(FILE *file)
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		const ParticleStore& atomparam = wspace.atom;
		SnapShotAtom *atom = wspace.cur.atom;

		int i;

		fprintf(file, "%8d !NBOND: bonds\n", bond.size());

		for(i = 0; i < bond.size(); i++) {
			fprintf(file, " %7d %7d", bond[i].i + 1, bond[i].j + 1);
			if((i % 4) == 3)
				fprintf(file, "\n");

		}

		fprintf(file, "\n\n");
		fprintf(file, "%8d !NTHETA: angles\n", angle.size());

		for(i = 0; i < angle.size(); i++) {
			fprintf(file, " %7d %7d %7d", angle[i].i + 1, angle[i].a + 1, angle[i].j + 1);
			if((i % 3) == 2)
				fprintf(file, "\n");
		}

		fprintf(file, "\n\n");
		fprintf(file, "%8d !NPHI: dihedrals\n", torsion.size());

		for(i = 0; i < torsion.size(); i++) {
			fprintf(file, " %7d %7d %7d %7d", torsion[i].i + 1, torsion[i].a + 1, torsion[i].b + 1, torsion[i].j + 1);
			if((i % 2) == 1)
				fprintf(file, "\n");
		}

		fprintf(file, "\n\n");

		fprintf(file, "%8d !NIMPHI: impropers\n", improper.size());

		for(i = 0; i < improper.size(); i++) {
			fprintf(file, " %7d %7d %7d %7d", improper[i].i + 1, improper[i].a + 1, improper[i].b + 1, improper[i].j + 1);
			if((i % 2) == 1)
				fprintf(file, "\n");
		}

		fprintf(file, "\n\n");

		fprintf(file, "%8d !NDON: donors\n", 0);
		fprintf(file, "\n\n");

		fprintf(file, "%8d !NACC: acceptors\n", 0);
		fprintf(file, "\n\n");

		fprintf(file, "%8d !NNB\n\n", 0);

		for(i = 0; i < wspace.atom.size(); i++) {
			fprintf(file, "%8d", 0);
			if((i % 8) == 7)
				fprintf(file, "\n");
		}
		fprintf(file, "\n\n");
		fprintf(file, "%8d %7d !NGRP\n%8d%8d%8d\n", 1, 0, 0, 0, 0);
		fprintf(file, "\n");

		return 0;
	}


} // namespace Physics


