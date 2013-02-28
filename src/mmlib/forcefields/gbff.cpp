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

// Pre-compiled Header
#include "global.h"

// Self header includes
// class definitions for generalised born solvent model implementation
#include "forcefields/gbff.h" 

#include "forcefields/ffparam.h"
#include "workspace/workspace.h"
#include "workspace/neighbourlist.h"
#include "workspace/bondorder.h"
// namespace includes
using namespace Maths;

namespace Physics
{
	FF_GeneralizedBorn_Still::FF_GeneralizedBorn_Still( WorkSpace &newwspace )
			: FF_GeneralizedBorn( newwspace )
	{
		FastMode = true;
		DielectricSolvent = 80.0;
		DielectricSolute = 1.0;
		GbsaStillCutoff = 10;
		ExpApproxThreshold = -4.0;
		DielectricOffset = -0.090;
		BornRadiusOffset = 0.0;
		
		ForceSwitch = false;
		EnergySwitch = true;

		P1 = 0.073; // these are the original still et al parameters
		P2 = 0.921; // & may be overwritten by the parameterfile
		P3 = 6.211; // From (9) Still, W. C.; Tempczyk, A.; Hawley,
		P4 = 15.236; // R. C.; Hendrickson, T. J. Am. Chem. Soc. 1990, 112, 6127.
		P5 = 1.254;

		epot_pol = 0;
	}

	void FF_GeneralizedBorn_Still::setup()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		Dielectric = DielectricSolute;
		FF_NonBonded::setup();

		bornpremul = -PhysicsConst::halfeconv_joule * (1 / DielectricSolute - 1 / DielectricSolvent);

		GBtype.clear();
		GB_AtomType t1;
		t1.radius = -1;
		GBtype.resize(wspace.ffps().AtomType.size(),t1);

		GB_atom_param_still.clear();
		GB_Atom_Param_Still t2;
		t2.terms123 = -1.0;
		t2.Vj = -1.0;
		GB_atom_param_still.resize( wspace.atom.size(), t2 );
		if(readGeneralisedBornSolvationSection() != 0)
		{
			throw(ProcedureException("Error reading Solvation parameters")); 
		}

		// prepare fast structure for lightning speed gbsa :)
		// it contains all the parameters in an AoS structure including
		// the non bonded parameters for vdw and vacuo elecstat

		GB_atom_param.clear();
		GB_Atom_Param t3;
		for(int i = 0; i < wspace.atom.size(); i++) 
		{
			t3.radiusij = wspace.atom[i].radius;
			t3.charge = wspace.atom[i].charge;
			t3.epsilon = wspace.atom[i].epsilon;
			GB_atom_param.push_back( t3 );
		}

		wspace.nlist().requestCutoff(Cutoff);

		calcFixedBornRadiiTerms();
		calcBornRadii_PairwiseApprox();

		wspace.nlist().calcShadow(true); // need shadow neighbors !!

		Active = true;
	}

	void FF_GeneralizedBorn_Still::info() const
	{
		// prints a little block of parameter information
		printf(" --Generalised Born Forcefield------- \n");
		if(!Active) printf(" -----------------------<INACTIVE>---\n");
		FF_NonBonded::info();
		printf(" Fast Mode:          %s\n", FastMode ? "yes" : "no");
		printf(" Born Radii:         Still et al.\n");
		printf(" Solvent Dielectric: %lf \n", DielectricSolvent);
		printf(" Solute Dielectric:  %lf \n", DielectricSolute);
		printf(" Still Cutoff:       %lf \n", GbsaStillCutoff);
		printf(" Dielec. Offset:     %lf \n", DielectricOffset);
		printf(" Exp Threshold:      %lf \n", ExpApproxThreshold);
	}

	void FF_GeneralizedBorn_Still::infoLine() const 
	{
		// prints a line of current energies
		const WorkSpace& wspace = getWSpace();
		const Hamiltonian *printene = &wspace.ene;
		if(FastMode) {
			printf("\t%6.1lf\t%6.1lf",
				double (printene->epot_vdw) * PhysicsConst::J2kcal * PhysicsConst::Na, double (printene->epot_elec) * PhysicsConst::J2kcal * PhysicsConst::Na);
		}
		printf("\t%6.1lf\t%6.1lf\t%6.1lf\t%7.2lf",
			epot_pol_self * PhysicsConst::J2kcal * PhysicsConst::Na,
			epot_pol_cross * PhysicsConst::J2kcal * PhysicsConst::Na, epot_pol * PhysicsConst::J2kcal * PhysicsConst::Na, (epot_pol + double (wspace.ene.epot_elec))
			*PhysicsConst::J2kcal * PhysicsConst::Na);

	}

	void FF_GeneralizedBorn_Still::infoLineHeader() const 
	{
		// prints the headers for the above function
		if(FastMode) {
			printf("\t%6s\t%6s", "EVdw", "Ecoul");
		}
		printf("\t%6s\t%6s\t%6s\t%7s", "Eself", "Ecross", "Epol", "Epol+ES");
	}

	void FF_GeneralizedBorn_Still::calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level)
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		if(level == Summary) {
			calcForces();
			printf(" VdW:           %10.3lf kcal/mol\n", double (wspace.ene.epot_vdw) * PhysicsConst::J2kcal * PhysicsConst::Na);
			printf(" Electrostatic: %10.3lf kcal/mol\n", double (wspace.ene.epot_elec) * PhysicsConst::J2kcal * PhysicsConst::Na);
			printf(" GB:            %10.3lf kcal/mol\n", epot_pol * PhysicsConst::J2kcal * PhysicsConst::Na);
			return;
		}

		calcBornRadii_PairwiseApprox();
		FF_NonBonded::calcEnergiesVerbose(level);
		calcBornEnergy_verbose(level);
		calcBornEnergy();

		const NeighbourData *fnbor = wspace.nlist().getData();
		printf(" Nr Neigbrs alpha Vj Term123 Charge self offdiag epol \n\n");
		for(int i = 0; i < wspace.atom.size(); i++) {
			printf("GB_atom_param_still: %4d %6d | %6.3lf %7.3lf %7.3lf %6.3lf \n",
				i,
				fnbor[i].n,
				GB_atom_param[i].bornradius,
				GB_atom_param_still[i].Vj,
				GB_atom_param_still[i].terms123,
				wspace.atom[i].charge);
		}
	}


	void FF_GeneralizedBorn_Still::calcEnergies()
	{
		calcBornRadii_PairwiseApprox();
		FF_NonBonded::calcEnergies();
		calcBornEnergy();
	}

	void FF_GeneralizedBorn_Still::calcForces()
	{
		calcBornRadii_PairwiseApprox();
		if(FastMode) 
		{
			calcForcesIncludingVacuo();
		} 
		else 
		{
			FF_NonBonded::calcForces();
			calcBornForces();
		}
	}


	int FF_GeneralizedBorn_Still::readGeneralisedBornSolvationSection()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		int line, attype;
		int errorstatus = 0;
		unsigned i;
		std::string sectionname = "SOLVATION_GENERALISED_BORN";
		Section section;
		if(wspace.ffps().getSection(sectionname,section)!=0){
			printf("ERROR: Cannot find section %s in forcefield parameters\n",
				sectionname.c_str());
			return -1;
		}

		int radiiGiven = true;
		int generalGiven = false; // have *all* specific parameters been given ?

		double myP1 = -1, myP2 = -1, myP3 = -1, myP4 = -1, myP5 = -1;// new general parameters
		// double newP1,newP2,newP3,newP4,newP5; // dummy Temporary variables

		//Start reading parameter file

		line = section.lineoffset;
		for(i=0;i<section.sectionline.size();i++){
			std::string linestring = section.sectionline[i];
			removecomments(linestring,"#\12\15");
			line++;
			std::vector<std::string> token;
			token = chopstr(linestring," \12\15\t");
			if(token.size() <= 0) continue; // empty line
			std::string command = token[0];


			if(cmpstring(command, "GENERALPARAMETERS")) { // TYPE syntax reading/checking of general parameters
				if(token.size() < 6) {
					printf("SYNTAX ERROR (line %d): Insufficient parameters for GENERALPARAMETERS\n", line);
					errorstatus = 1;
					continue;
				}

				if(str2double(token[1],P1)!=0){ // assign parameters
					printf("ERROR: Numerical value expected for parameter P1 after GENERALPARAMETERS\n");
					errorstatus = 1;
					continue;
				}
				if(str2double(token[2],P2)!=0){ // assign parameters
					printf("ERROR: Numerical value expected for parameter P2 after GENERALPARAMETERS\n");
					errorstatus = 1;
					continue;
				}
				if(str2double(token[3],P3)!=0){ // assign parameters
					printf("ERROR: Numerical value expected for parameter P3 after GENERALPARAMETERS\n");
					errorstatus = 1;
					continue;
				}
				if(str2double(token[4],P4)!=0){ // assign parameters
					printf("ERROR: Numerical value expected for parameter P4 after GENERALPARAMETERS\n");
					errorstatus = 1;
					continue;
				}
				if(str2double(token[5],P5)!=0){ // assign parameters
					printf("ERROR: Numerical value expected for parameter P5 after GENERALPARAMETERS\n");
					errorstatus = 1;
					continue;
				}

				generalGiven = true; // take note that the general parameters were implicitly defined
			} else
				if(cmpstring(command, "SOLVATIONTYPE")) { // TYPE syntax reading/checking of SOLVATIONTYPE
					if(token.size() < 3) {
						printf
							("SYNTAX ERROR (line %d): Insufficient parameters for SOLVATIONTYPE (either 2 or 7 parameters required)\n",
							line);
						errorstatus = 1;
						continue;
					}
					attype = wspace.ffps().findAtomType(token[1]); // find attype number
					if(attype < 0) { // in case we didn't find it
						printf("SYNTAX ERROR (line %d): Atom name %s unknown \n", line, token[1].c_str());
						errorstatus = 1;
						continue;
					}

					if(str2double(token[2],GBtype[attype].radius)!=0){ // assign parameters
						printf("ERROR: Numerical value expected for parameter 2 after SOLVATIONTYPE\n");
						errorstatus = 1;
						continue;
					}

				} else {
					printf("ERROR: Identifier %s not known in %s section \n", command.c_str(), sectionname.c_str());
					errorstatus = 1;
				}
		}


		// Error handling ----------------------------------------------------

		// test if *all* of the radii were given, otherwise falsifiy specificsGivenFlag
		for(i = 0; i < wspace.ffps().AtomType.size(); i++) {
			if(GBtype[i].radius < 0){
				fprintf(stderr,"ERROR: No Parameters found for atom Type %s \n",wspace.ffps().AtomType[i].name.c_str());
				radiiGiven = false;
			}
		}
		if(!radiiGiven) {
			printf("ERROR: Not all GBSA Radii were found in forcefield file !\n");
		}
		// 'handle' error status
		if(errorstatus != 0) THROW(ProcedureException,"Errors during readout of Generalised Born solvation parameters");

		printf("ff.gbsa: Using FF_GeneralizedBorn_Still parameters: %lf %lf %lf %lf %lf ", P1, P2, P3, P4, P5);
		if(generalGiven)			printf(" - as specifed\n");
		else									printf(" - by default, not implicitly defined (warning)\n");

		this->needsetup = false;
		return errorstatus;
	}



	// calculates the fixed terms (depended on connectivity not depended on coordinates)
	int FF_GeneralizedBorn_Still::calcFixedBornRadiiTerms()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		int i, nj, j, nk, k;
		double Rj, Rk, hjk, rjk;
		double rij;

		double term1, term2, term3;
		int Type;

		wspace.nlist().calcNewList();

		for(j = 0; j < wspace.atom.size(); j++) {
			Type = wspace.atom[j].FFType;
			Rj = GBtype[Type].radius; //wspace.atom[i].radius;

			GB_atom_param_still[j].Vj = 0;
			if(Rj <= 0.01)
				Rj = 1.15; // according to Qui et al. for hydrogen atoms

			GB_atom_param_still[j].Vj = Maths::MathConst::FourPI * Rj * Rj * Rj / 3; // its own volume

		  //printf("Rad %5d %5.2lf Volume %8.3lf ",j,Rj,GB_atom_param_still[j].Vj);
			for(nk = 0; nk < wspace.atom[j].cov12atom.size(); nk++) {// minus half the overlapping volumes of
				// all 1,2 bonded neighbours
				k = wspace.atom[j].cov12atom[nk].i;

				Rk = GBtype[wspace.atom[k].FFType].radius;
				if(Rk <= 0.01)
					Rk = 1.15; // according to Qui et al. for hydrogen atoms
				rjk = 1.01 * wspace.atom[j].cov12atom[nk].d;
				hjk = Rj * (1.0 + (sqr(Rk) - sqr(Rj) - sqr(rjk)) / (2.0 * Rj * rjk));

				GB_atom_param_still[j].Vj -= Maths::MathConst::PI * sqr(hjk) * (3 * Rj - hjk) / 3;
			}

			if(GB_atom_param_still[j].Vj < 0)
				GB_atom_param_still[j].Vj = 0;
		}

		for(i = 0; i < wspace.atom.size(); i++) {
			GB_atom_param_still[i].terms123 = 0;
			Type = wspace.atom[i].FFType;
			Rj = GBtype[Type].radius; 

			if(Rj <= 0.01)
				Rj = 1.15; // according to Qui et al. for hydrogen atoms

			// see [2] for methods
			//constant self terms
			term1 = -166.0 / (Rj + DielectricOffset + P1);
			//constant 1-2 atom terms
			term2 = 0;
			for(nj = 0; nj < wspace.atom[i].cov12atom.size(); nj++) {
				j = wspace.atom[i].cov12atom[nj].i;
				rij = wspace.atom[i].cov12atom[nj].d;
				
				term2 += P2 * GB_atom_param_still[j].Vj / (sqr(rij) * sqr(rij));
			}
			//constant 1-3 atom terms
			term3 = 0;
			for(nj = 0; nj < wspace.atom[i].cov13atom.size(); nj++) {
				j = wspace.atom[i].cov13atom[nj].i;
				rij = wspace.atom[i].cov13atom[nj].d;
				term3 += P3 * GB_atom_param_still[j].Vj / (sqr(rij) * sqr(rij));
			}

			GB_atom_param_still[i].terms123 = term1 + term2 + term3;
			GB_atom_param[i].bornradius = -PhysicsConst::halfeconv / GB_atom_param_still[i].terms123;// born radius so far
		}

		return 0;
	}


	// -------------------------------------------------------------------------------
	// Class: FF_GeneralizedBorn_Still
	// Function: calcBornRadii_PairwiseApprox();
	// -------------------------------------------------------------------------------
	// Parameters: none
	//
	// This implements the pairwise approximation of bornradii using the Clarke Still
	// approximation. The derivatives of the bornradii are calculated in calcBornForces();
	// or calcBornForces_vacuo(); depending on mode

	int FF_GeneralizedBorn_Still::calcBornRadii_constant(double constBornRadius)
	{
		WorkSpace& wspace = getWSpace();
		for(int i = 0; i < wspace.atom.size(); i++) {
			if(GB_atom_param[i].bornradius<constBornRadius)
				GB_atom_param[i].bornradius = constBornRadius;
		}
		return 0;
	}


	// -------------------------------------------------------------------------------
	// Class: FF_GeneralizedBorn_Still
	// Function: calcBornRadii_PairwiseApprox();
	// -------------------------------------------------------------------------------
	// Parameters: none
	//
	// This implements the pairwise approximation of bornradii using the Clarke Still
	// approximation. The derivatives of the bornradii are calculated in calcBornForces();
	// or calcBornForces_vacuo(); depending on mode

	int FF_GeneralizedBorn_Still::calcBornRadii_PairwiseApprox()
	{
		WorkSpace& wspace = getWSpace();

		int i, nj, j;
		int Type;
		double Gpol;
		double Dist_ij;
		double invP5;
		double CCF;
		double rdivR;
		double Rj, Ri;

		const NeighbourData *fnbor = wspace.nlist().getData();

		invP5 = 1 / P5;

		for(i = 0; i < wspace.atom.size(); i++) {
			Type = wspace.atom[i].FFType;
			Ri = GBtype[Type].radius;
			// double count=0.0;

			Gpol = GB_atom_param_still[i].terms123;

			//run over all neighbors
			for(nj = 0; nj < fnbor[i].n; nj++) {
				j = NList32Bit_Index(fnbor[i].i[nj]);
				Dist_ij = sqrdist(wspace.cur.atom[i].p,wspace.cur.atom[j].p);

				if(Dist_ij > sqr(GbsaStillCutoff))continue;
			


				if( NList32Bit_BondOrder(fnbor[i].i[nj]) <= BondOrder_1_3_Pair)
					continue;// only get non-bonded and 1-4 neighbors
		

				Rj = GBtype[wspace.atom[j].FFType].radius; //wspace.atom[i].radius;

				rdivR = Dist_ij / sqr((Ri + Rj));
				if(rdivR > invP5)
					CCF = 1.0;
				else
					CCF = sqr(0.5 - 0.5 * cos((rdivR * P5 * Maths::MathConst::PI)));

				Gpol += P4 * GB_atom_param_still[j].Vj * CCF / (sqr(Dist_ij));
				// count += 1.0;
			}
			GB_atom_param[i].bornradius = -PhysicsConst::halfeconv / Gpol + BornRadiusOffset;
		}
		return 0;
	}




	void FF_GeneralizedBorn_Still::calcBornEnergy_covalentTerms()
	{
		WorkSpace& wspace = getWSpace();

		int i, nj, j, atoms;
		double aiaj, srij;
		double epol;
		double epot_elecstat;

		epot_covalent = 0;
		epot_elecstat = 0;

		atoms = wspace.atom.size();

		for(i = 0; i < atoms; i++) {
			for(nj = 0; nj < wspace.atom[i].cov12atom.size(); nj++) {
				j = wspace.atom[i].cov12atom[nj].i;

				epol = wspace.atom[i].charge * wspace.atom[j].charge;

				srij = sqr(wspace.atom[i].cov12atom[nj].d);

				aiaj = GB_atom_param[i].bornradius * GB_atom_param[j].bornradius;
				epol /= sqrt(srij + aiaj * exp(-srij * 0.25 / aiaj));

				epot_elecstat += epol;
			}

			for(nj = 0; nj < wspace.atom[i].cov13atom.size(); nj++) {
				j = wspace.atom[i].cov13atom[nj].i;

				epol = wspace.atom[i].charge * wspace.atom[j].charge;

				srij = sqr(wspace.atom[i].cov13atom[nj].d);

				aiaj = GB_atom_param[i].bornradius * GB_atom_param[j].bornradius;
				epol /= sqrt(srij + aiaj * exp(-srij * 0.25 / aiaj));

				epot_elecstat += epol;
			}
		}

		epot_elecstat *= bornpremul;
		epot_covalent = epot_elecstat;

		printf("epot_covalent = %lf \n", epot_covalent);
	}


	void FF_GeneralizedBorn_Still::calcBornEnergy()
	{
		WorkSpace& wspace = getWSpace();

		int i, nj, j, atoms;
		double ai, qi, qj, aiaj, srij;
		double epol;
		double Dist_ij;
		double Swidth = Cutoff - InnerCutoff;
		double invSwidth = 1 / Swidth;
		double S;

		const NeighbourData *fnbor = wspace.nlist().getData();
		atoms = wspace.atom.size();

		epot_pol = 0;
		epot_pol_self = 0.0;
		epot_pol_cross = 0.0;		

		for(i = 0; i < atoms; i++) {
			qi = wspace.atom[i].charge;
			ai = GB_atom_param[i].bornradius;
			epot_pol -= sqr(qi) / ai;
			epot_pol_self -= sqr(qi) / ai;

			// do neighbors
			for(nj = 0; nj < fnbor[i].n; nj++) {
				j = NList32Bit_Index(fnbor[i].i[nj]);
				if(i > j)
					continue; // only calculate half the matrix (fGB is symetric)

				Dist_ij = dist(wspace.cur.atom[i].p,wspace.cur.atom[j].p);
				if(Dist_ij>Cutoff) continue;

				// Switching function (S)
				if(Dist_ij < InnerCutoff) {
					S = 1.0;
				} else {
					S = (1.0 - sqr(invSwidth * (Dist_ij - InnerCutoff)));
					S = sqr(S);
				}

				srij = sqr(Dist_ij);

				qj = wspace.atom[j].charge;

				aiaj = ai * GB_atom_param[j].bornradius;
				epol = 2.0 * S * qi * qj / sqrt(srij + aiaj * exp(-srij * 0.25 / aiaj));

				epot_pol -= epol;
				epot_pol_cross -= epol;
			}
		}

		epot_pol *= -bornpremul;
		epot_pol_cross *= -bornpremul;
		epot_pol_self *= -bornpremul;

		wspace.ene.epot += epot_pol;
	}


	void FF_GeneralizedBorn_Still::calcBornEnergy_verbose(AtomicVerbosity verboselevel)
	{
		WorkSpace& wspace = getWSpace();

		int i, nj, j, atoms;
		double aiaj, srij;
		double epol;
		double epot_elecstat;
		double epot_elecvac, elecvac;

		double Dist_ij;
		double Swidth = Cutoff - InnerCutoff;
		double invSwidth = 1 / Swidth;
		double S;

		epot_pol = 0;
		epot_elecstat = 0;
		epot_elecvac = 0;

		const NeighbourData *fnbor = wspace.nlist().getData();
		atoms = wspace.atom.size();

		for(i = 0; i < atoms; i++) {
			for(nj = 0; nj < fnbor[i].n; nj++) {// do neighbors
				j = NList32Bit_Index(fnbor[i].i[nj]);

				Dist_ij = dist(wspace.cur.atom[i].p,wspace.cur.atom[j].p);

				// Switching function (S)
				if(Dist_ij < InnerCutoff) {
					S = 1.0;
				} else {
					S = (1.0 - sqr(invSwidth * (Dist_ij - InnerCutoff)));
					S = sqr(S);
				}

				srij = sqr(Dist_ij);

				epol = wspace.atom[i].charge * wspace.atom[j].charge;

				elecvac = epol / sqrt(srij);

				aiaj = GB_atom_param[i].bornradius * GB_atom_param[j].bornradius;
				epol /= sqrt(srij + aiaj * exp(-srij * 0.25 / aiaj));
				//epol /= sqrt(srij + aiaj);

				epol *= S; // apply switch function
				elecvac *= S;

				epot_elecstat += epol;

				if(verboselevel == Detailed) {
					printf("Solv:%5d(%4s)%5d(%4s) %5.2lf % 4.2lf % 4.2lf % 8.3lf % 8.3lf %6.3lf\n",
						i, wspace.atom[i].pdbname.c_str(),
						j, wspace.atom[j].pdbname.c_str(),
						sqrt(srij),
						wspace.atom[i].charge,
						wspace.atom[j].charge,
						bornpremul * epol, bornpremul * elecvac, elecvac / (elecvac - bornpremul * epol));
				}

			}
		}

		epot_elecstat *= bornpremul;

		for(i = 0; i < atoms; i++) {
			epol = sqr(wspace.atom[i].charge);
			epol /= GB_atom_param[i].bornradius;
			epot_pol += epol;
			epol *= bornpremul;
		}

		printf("Born Pre-Multiplication Factor: %10.6lf\n", bornpremul);

		epot_pol *= bornpremul;
		epot_pol += epot_elecstat;
	}










	// -------------------------------------------------------------------------------
	// Class: FF_GeneralizedBorn_Still
	// Function: calcBornForces_numerical(double dc, int dtype);
	// -------------------------------------------------------------------------------
	// Parameters: double dc Numerical Differentitation Displacement
	// int dtype Differential Type: 0 - full, all born radii change
	// 1 - only ai changes (aj is constant)
	// 2 - all born radii are constant
	// Purpose: This calculates numerical derivatives of the born energies, by
	// taking each atom in turn, moving it in all 3 space dimensions
	// and calculating the change in total polarisation energy with respect to
	// atomic displacements. 3 Types possible , Type 0 is the true, rigorous
	// differential, the others are approximations
	//
	void FF_GeneralizedBorn_Still::calcBornForces_numerical(double dc, int dtype)
	{
		WorkSpace& wspace = getWSpace();

		int i, atoms, k, u;

		double epot_elecstat, epot_elecvac;

		double epol_0, epol_dx, epol_dy, epol_dz;

		dvector atomsav;

		dvector vdx(dc, 0.0, 0.0); // the dvector  displacements
		dvector vdy(0.0, dc, 0.0);
		dvector vdz(0.0, 0.0, dc);
		dvector dGpoldc;

		dvector force;


		epot_pol = 0;
		epot_elecstat = 0;
		epot_elecvac = 0;

		atoms = wspace.atom.size();

		double *a0 = new double[atoms];
		dvector *dadx = new dvector[atoms];

		calcBornRadii_PairwiseApprox(); // calcualte the born polarisation energy
		for(i = 0; i < atoms; i++)
			a0[i] = GB_atom_param[i].bornradius;

		calcBornEnergy(); // for the initial position of all atoms
		epol_0 = epot_pol; // and save it in epol_0

		for(i = 0; i < atoms; i++) {

			atomsav.setTo(wspace.cur.atom[i].p);
			for(k = 0; k < 3; k++) {
				// reset moved atom to its original position
				switch (k) 
				{
				case 0:
					wspace.cur.atom[i].p.add(vdx); // add a small increment
					calcBornRadii_PairwiseApprox(); // calculate the born polarisation energy

					for(u = 0; u < atoms; u++)
						dadx[u].x = (a0[u] - GB_atom_param[u].bornradius) / dc;

					if(dtype > 0) {
						for(u = 0; u < atoms; u++) {
							if(dtype <= 1)
								if(u == i)
									continue;
							GB_atom_param[u].bornradius = a0[u];
						}
					}

					calcBornEnergy(); // for the new position of all atoms
					epol_dx = epot_pol - epol_0;
					break;
				case 1:
					wspace.cur.atom[i].p.add(vdy); // add a small increment
					calcBornRadii_PairwiseApprox(); // calcualte the born polarisation energy

					for(u = 0; u < atoms; u++)
						dadx[u].y = (a0[u] - GB_atom_param[u].bornradius) / dc;

					if(dtype > 0) {
						for(u = 0; u < atoms; u++) {
							if(dtype <= 1)
								if(u == i)
									continue;
							GB_atom_param[u].bornradius = a0[u];
						}
					}

					calcBornEnergy(); // for the new position of all atoms
					epol_dy = epot_pol - epol_0;
					break;
				case 2:
					wspace.cur.atom[i].p.add(vdz); // add a small increment
					calcBornRadii_PairwiseApprox(); // calcualte the born polarisation energy

					for(u = 0; u < atoms; u++)
						dadx[u].z = (a0[u] - GB_atom_param[u].bornradius) / dc;

					if(dtype > 0) {
						for(u = 0; u < atoms; u++) {
							if(dtype <= 1)
								if(u == i)
									continue;
							GB_atom_param[u].bornradius = a0[u];
						}
					}

					calcBornEnergy(); // for the new position of all atoms
					epol_dz = epot_pol - epol_0;
					break;
				}
				wspace.cur.atom[i].p.setTo(atomsav);
			}

			// for(u=0;u<atoms;u++) printf("\t%5.3lf",dadx[u].mag());
			// printf("\n");
			// printf("\t%8.6lf\n",dadx[i].mag());
			//

			force.setTo(epol_dx / dc, epol_dy / dc, epol_dz / dc);
			force.mul(-1 / (PhysicsConst::Angstrom));
			wspace.cur.atom[i].f.add(force);
		}

		for(i = 0; i < atoms; i++) {
			wspace.old.atom[i].f.setTo(wspace.cur.atom[i].f);
			wspace.cur.atom[i].f.zero();
		}

		// just redo the clean calculation to finish it off
		calcBornRadii_PairwiseApprox(); // calcualte the born polarisation energy
		calcBornEnergy(); // for the initial position of all atoms

		wspace.ene.epot += epot_pol;

		delete[]a0;
		delete[]dadx;
	}







	// -------------------------------------------------------------------------------
	// Class: FF_GeneralizedBorn_Still
	// Function: calcBornForces_constantBornRadii();
	// -------------------------------------------------------------------------------
	// Parameters: none
	//
	// This is the simplest implementation of the Bornradii derivatives
	// it (foolishly) assumes that the bornradii remain constant on infinitesimal
	// movements of the atoms - this makes the derivation considerably more managable
	// in terms of speed and simplicity - yet looses out on accuracy.
	// especially the selfenergies will feel no force, since they are solely dependent
	// on the bornradii not on interactomic distances - hence assuming the former
	// are constants, the derivatives off the self terms will be 0 (and hence not calculated)
	// Yet, the total energy contribution of any given atom depends largely on the
	// cross terms, and thus the derivative optained in this fashion vary only
	// slightly from their true values (see comparison). This derivative corresponds to
	// the Type 1 derivative of the numerical implementation

	void FF_GeneralizedBorn_Still::calcBornForces_constantBornRadii()
	{
		int i, nj, j, atoms;
		double aiaj, Dist_ij, srij;
		double sqrt_lGB, lGB;
		double qiqj, epol;
		double epot_elecstat;
		double deijdrij2;

		double Swidth = Cutoff - InnerCutoff;
		double invSwidth = 1 / Swidth;
		double S, dSdrij2;

		epot_pol = 0;
		epot_pol_self = 0.0;
		epot_pol_cross = 0.0;
		epot_elecstat = 0;

		// Proxies
		WorkSpace& wspace = getWSpace();
		const NeighbourData *fnbor = wspace.nlist().getData();
		atoms = wspace.atom.size();

		for(i = 0; i < atoms; i++) {
			GB_atom_param[i].dedi.setTo(0, 0, 0);

			for(nj = 0; nj < fnbor[i].n; nj++) {// do neighbors
				j = NList32Bit_Index(fnbor[i].i[nj]);
				Dist_ij = dist(wspace.cur.atom[i].p,wspace.cur.atom[j].p);

				// Switching function (S)
				if(Dist_ij < InnerCutoff) {
					S = 1;
					dSdrij2 = 0;
				} else {
					S = (1 - sqr(invSwidth * (Dist_ij - InnerCutoff)));
					dSdrij2 = S * 2 * sqr(invSwidth) * (Dist_ij - InnerCutoff) / Dist_ij;
					S = sqr(S);

				}

				qiqj = wspace.atom[i].charge * wspace.atom[j].charge;

				srij = sqr(Dist_ij);

				aiaj = GB_atom_param[i].bornradius * GB_atom_param[j].bornradius;

				lGB = srij + aiaj * exp(-srij * 0.25 / aiaj);
				sqrt_lGB = sqrt(lGB);
				epol = qiqj / sqrt_lGB;

				deijdrij2 = -0.5 * qiqj / (lGB * sqrt_lGB) * (1 - 0.25 * exp(-srij * 0.25 / aiaj));
				deijdrij2 = S * deijdrij2 + epol * dSdrij2;
				epol *= S;

				GB_atom_param[i].dedi.x += deijdrij2 * 2 * (wspace.cur.atom[i].p.x - wspace.cur.atom[j].p.x);
				GB_atom_param[i].dedi.y += deijdrij2 * 2 * (wspace.cur.atom[i].p.y - wspace.cur.atom[j].p.y);
				GB_atom_param[i].dedi.z += deijdrij2 * 2 * (wspace.cur.atom[i].p.z - wspace.cur.atom[j].p.z);



				epot_elecstat += epol;
				epot_pol_cross += epol;
			}
			GB_atom_param[i].dedi.mul(2 * bornpremul);

		}

		epot_elecstat *= bornpremul;

		// for self terms only calculate the energies, the derivatives will be
		// 0 anyway
		dvector force;
		for(i = 0; i < atoms; i++) {
			force.setTo(GB_atom_param[i].dedi);
			force.mul(-1 / (PhysicsConst::Angstrom));
			wspace.cur.atom[i].f.add(force);

			epol = sqr(wspace.atom[i].charge);
			epol /= GB_atom_param[i].bornradius;
			epot_pol += epol;
			epol *= bornpremul;
		}

		epot_pol *= bornpremul;
		epot_pol_cross *= bornpremul;
		epot_pol_self *= bornpremul;
		epot_pol += epot_elecstat;
	}




	void FF_GeneralizedBorn_Still::calcForces_constantBornRadii_IncludingVacuo()
	{
		// set up proxies to workspace to make code more readable.
		WorkSpace& wspace = getWSpace();
		size_t   atoms = wspace.atom.size();                 // number of atoms in workspace
		SnapShotAtom *atom = wspace.cur.atom; // atom coordinate array

		int i, j, nj;

		int typej;

		double dedr;
		double dedx, dedy, dedz;

		double dddx, dddy, dddz;

		const NeighbourData *fnbor = wspace.nlist().getData();

		double Swidth = Cutoff - InnerCutoff;
		double invSwidth = 1 / Swidth;
		double S;
		double dSdd;

		double Dist_ij, invdistij;
		double qi, ai, aj, aiaj, srij;
		double uGB, fGB, lGB, mGB, tGB;
		double premul;
		double qiqj;
		double pair_elec14scaling = 1.0;

		double epol;
		double dedai;

		double atomix, atomiy, atomiz, force_magnitude;
		dvector fv, force;
		double vdw_force = 0, // individual force magnitudes
			elec_force = 0;
		double vdw_potential = 0, // individual potential energy contributions
			elec_potential = 0;

		double pair_vdw14scaling = 1.0;
		double radiusij, epsilon, A, B;
		double invdielectric = 1 / DielectricSolute;
		double invsolventdielectric = 1 / DielectricSolvent;

		double atomi_radius;
		double atomi_epsilon;
		dvector atomi_p;

		wspace.ene.epot_vdw = 0;
		wspace.ene.epot_elec = 0;

		epot_pol = 0.0;
		epot_pol_self = 0.0;
		epot_pol_cross = 0.0;

		atoms = wspace.atom.size();

		for(i = 0; i < atoms; i++)
			GB_atom_param[i].position.setTo(atom[i].p);
		premul = bornpremul;

		for(i = 0; i < atoms; i++) {
			GB_atom_param[i].dedi.setTo(0, 0, 0);
		}
		for(i = 0; i < atoms; i++) {
			ai = GB_atom_param[i].bornradius;
			atomix = atom[i].p.x;
			atomiy = atom[i].p.y;
			atomiz = atom[i].p.z;
			dedai = 0;

			atomi_radius = GB_atom_param[i].radiusij;
			atomi_epsilon = GB_atom_param[i].epsilon;
			qi = GB_atom_param[i].charge;

			// do neighbors
			for(nj = 0; nj < fnbor[i].n; nj++) {
				j = NList32Bit_Index(fnbor[i].i[nj]);
				if(i > j)
					break;

				typej = NList32Bit_BondOrder(fnbor[i].i[nj]);


				aj = GB_atom_param[j].bornradius;
				qiqj = qi * GB_atom_param[j].charge;

				dddx = (double)(atomix - GB_atom_param[j].position.x); // these are incomplete (need div by Dist_ij)
				dddy = (double)(atomiy - GB_atom_param[j].position.y);
				dddz = (double)(atomiz - GB_atom_param[j].position.z);
				srij = sqr(dddx) + sqr(dddy) + sqr(dddz);
				Dist_ij = sqrt(srij);

				force_magnitude = 0; // total force/total pot. energy = 0

				invdistij = 1 / Dist_ij; // inverse distance - universally required

				aiaj = ai * aj;

				uGB = -srij * 0.25 / aiaj;

				pair_elec14scaling = 1.0;
				if(NList32Bit_BondOrder(fnbor[i].i[nj]) == BondOrder_1_4_Pair)
					pair_elec14scaling = Elec14Scaling;

				// exact GB and vacuo electrostatic
				if(uGB > ExpApproxThreshold) {
					tGB = exp(uGB);
					mGB = tGB * aiaj;
					lGB = srij + mGB;

					fGB = sqrt(lGB);
					epol = premul * qiqj / fGB;

					dedr = -epol * (Dist_ij - 0.25 * Dist_ij * tGB) / lGB;

					elec_potential = 0;
					elec_force = 0;
					if(typej <= 1) { // only do vdw/elec for non-bonded and 1-4 neighbors
						elec_potential = (double) PhysicsConst::econv_joule *invdielectric * pair_elec14scaling * qiqj * invdistij;
						elec_force = -elec_potential * invdistij;
					}
					// Apply switching function (S)
					if(Dist_ij >= InnerCutoff) {
						S = (1.0 - sqr(invSwidth * (Dist_ij - InnerCutoff))); //Dimension less
						dSdd = -4.0 * sqr(invSwidth) * (Dist_ij - InnerCutoff) * S; //Units of A-1
						S = sqr(S);

						elec_force = S * elec_force + elec_potential * dSdd;
						elec_potential *= S;
						dedr = S * dedr + dSdd * epol;
						epol *= S;
					}

					dedr *= 2 * invdistij;
					dedx = dedr * dddx;
					dedy = dedr * dddy;
					dedz = dedr * dddz;

					epot_pol += epol;
					epot_pol += epol;
					epot_pol_cross += epol;
					epot_pol_cross += epol;

					GB_atom_param[i].dedi.x += dedx;
					GB_atom_param[i].dedi.y += dedy;
					GB_atom_param[i].dedi.z += dedz;
					GB_atom_param[j].dedi.x -= dedx;
					GB_atom_param[j].dedi.y -= dedy;
					GB_atom_param[j].dedi.z -= dedz;

					// water Dielectric electrostatics (long distance)
				} else {
					if(typej <= 1) { // only do vdw/elec for non-bonded and 1-4 neighbors
						elec_potential = PhysicsConst::econv_joule * invsolventdielectric * pair_elec14scaling * qiqj * invdistij;
						elec_force = -elec_potential * invdistij;

						// Apply switching function (S)
						if(Dist_ij >= InnerCutoff) {
							S = (1.0 - sqr(invSwidth * (Dist_ij - InnerCutoff))); //Dimension less
							dSdd = -4.0 * sqr(invSwidth) * (Dist_ij - InnerCutoff) * S; //Units of A-1
							S = sqr(S);
							elec_force = S * elec_force + elec_potential * dSdd;
							elec_potential *= S;
						}
					}
				}

				pair_vdw14scaling = 1.0;
				if(typej <= 1) { // only do vdw/elec for non-bonded and 1-4 neighbors
					if(typej == 1)
						pair_vdw14scaling = Vdw14Scaling;

					// Van der Waals Force
					radiusij = sqr(invdistij * (atomi_radius + GB_atom_param[j].radiusij));
					epsilon = pair_vdw14scaling * atomi_epsilon * GB_atom_param[j].epsilon;

					//invd6 = cube()*cube(invdistij); // 1/d^6 needed for vdw
					B = cube(radiusij);
					A = sqr(B);

					vdw_potential = (2.0 * epsilon) * (0.5 * A - B);
					vdw_force = -(12.0 * epsilon * invdistij) * (A - B);

					// add up potentials -----------------------------

					wspace.ene.epot += vdw_potential;
					wspace.ene.epot_vdw += vdw_potential;

					wspace.ene.epot += elec_potential;
					wspace.ene.epot_elec += elec_potential;

					// now apply the force
					force_magnitude= 0;
					force_magnitude+= vdw_force;
					force_magnitude+= elec_force;
					force_magnitude*= invdistij;
					dedx = force_magnitude* dddx;
					dedy = force_magnitude* dddy;
					dedz = force_magnitude* dddz;

					GB_atom_param[i].dedi.x += dedx;
					GB_atom_param[i].dedi.y += dedy;
					GB_atom_param[i].dedi.z += dedz;
					GB_atom_param[j].dedi.x -= dedx;
					GB_atom_param[j].dedi.y -= dedy;
					GB_atom_param[j].dedi.z -= dedz;
				}
			}
			epol = premul * sqr(qi) / ai;
			epot_pol += epol; // add up energies
			epot_pol_self += epol;
		}

		//add electrostatic energy to the total
		wspace.ene.epot += epot_pol;
		//add electrostatic energy to the energy components
		wspace.ene.epot_pol += epot_pol;
		wspace.ene.epot_pol_cross += epot_pol_cross;
		wspace.ene.epot_pol_self += epot_pol_self;

		for(i = 0; i < atoms; i++) {
			force.setTo(GB_atom_param[i].dedi);
			force.mul((double) - 1.0 / (double) (PhysicsConst::Angstrom));
			wspace.cur.atom[i].f.add(force);
		}
	}





	// -------------------------------------------------------------------------------
	// Class: FF_GeneralizedBorn_Still
	// Function: calcBornForces();
	// -------------------------------------------------------------------------------
	// Parameters: none
	//
	// This implements the true derivative using the born raddi chain rule of
	// partial derivatives with respect to rij, ai and all aj

	void FF_GeneralizedBorn_Still::calcBornForces(){
		int i, j, nj;
		int atoms;
		double dedr, dedalpha;
		double dedx, dedy, dedz;

		double dddx, dddy, dddz;

		double sqrcutoff = sqr(Cutoff);
		double Swidth = Cutoff - InnerCutoff;
		double invSwidth = 1 / Swidth;
		double S;
		double dSdd;

		double Dist_ij, invdistij;
		double qi, ai, aj, aiaj, srij;
		double uGB, fGB, lGB, mGB, tGB;
		double premul;
		double qiqj;

		double epol;
		double dedai;

		double atomix, atomiy, atomiz;
		dvector force;

		epot_pol = 0.0;
		epot_pol_self = 0.0;
		epot_pol_cross = 0.0;

		// Proxies
		WorkSpace& wspace = getWSpace();
		const NeighbourData *fnbor = wspace.nlist().getData();
		atoms = wspace.atom.size();

		premul = bornpremul;

		// set all arrays to 0
		for(i = 0; i < atoms; i++) {
			GB_atom_param[i].dedi.setTo(0, 0, 0);
			GB_atom_param[i].deda = 0.0;
		}

		int paircount=0;

		for(i = 0; i < atoms; i++) {
			ai = GB_atom_param[i].bornradius;
			qi = premul * wspace.atom[i].charge;

			atomix = wspace.cur.atom[i].p.x;
			atomiy = wspace.cur.atom[i].p.y;
			atomiz = wspace.cur.atom[i].p.z;
			dedai = 0;
			// do neighbors
			for(nj = 0; nj < fnbor[i].n; nj++) {

				j = NList32Bit_Index(fnbor[i].i[nj]);
				if(j > i) break; // ignore shadow neighbors

				//Dist_ij = fnbor[i].d[nj];
				//if(Dist_ij>Cutoff){
				// printf("%lf \n",Dist_ij);
				// continue;
				//}
				aj = GB_atom_param[j].bornradius;
				qiqj = qi * wspace.atom[j].charge;
				paircount++;
				dddx = (atomix - wspace.cur.atom[j].p.x); // these are incomplete (need div by Dist_ij)
				dddy = (atomiy - wspace.cur.atom[j].p.y);
				dddz = (atomiz - wspace.cur.atom[j].p.z);
				srij = sqr(dddx) + sqr(dddy) + sqr(dddz);
				if(srij>sqrcutoff) continue;
				Dist_ij = sqrt(srij);

				invdistij = 1 / Dist_ij;

				aiaj = ai * aj;

				uGB = -srij * 0.25 / aiaj;
				if(uGB > ExpApproxThreshold) {
					tGB = exp(uGB);
					mGB = tGB * aiaj;
					lGB = srij + mGB;

					fGB = sqrt(lGB);
					epol = qiqj / fGB;

					dedr = -epol * (Dist_ij - 0.25 * Dist_ij * tGB) / lGB;
					dedalpha = -epol * tGB * (1.0 + 0.25 * srij / aiaj) / lGB;

					// Apply switching function (S)
					if(Dist_ij >= InnerCutoff) {
						S = (1.0 - sqr(invSwidth * (Dist_ij - InnerCutoff))); //Dimension less
						dSdd = -4.0 * sqr(invSwidth) * (Dist_ij - InnerCutoff) * S; //Units of A-1
						S = sqr(S);

						dedr = S * dedr + dSdd * epol;
						dedalpha *= S;
						epol *= S;
					}
					dedr *= 2 * invdistij;
					dedx = dedr * dddx;
					dedy = dedr * dddy;
					dedz = dedr * dddz;

					epot_pol += epol;
					epot_pol += epol;
					epot_pol_cross += epol;
					epot_pol_cross += epol;

					GB_atom_param[i].dedi.x += dedx;
					GB_atom_param[i].dedi.y += dedy;
					GB_atom_param[i].dedi.z += dedz;
					GB_atom_param[j].dedi.x -= dedx;
					GB_atom_param[j].dedi.y -= dedy;
					GB_atom_param[j].dedi.z -= dedz;

					dedai += dedalpha * aj;
					GB_atom_param[j].deda += dedalpha * ai;
				} else {
					tGB = 0;
					mGB = 0;
					lGB = srij;
					fGB = Dist_ij;

					epol = qiqj * invdistij;
					dedr = -epol * invdistij;

					// Apply switching function (S)
					if(Dist_ij >= InnerCutoff) {
						S = (1.0 - sqr(invSwidth * (Dist_ij - InnerCutoff))); //Dimension less
						dSdd = -4.0 * sqr(invSwidth) * (Dist_ij - InnerCutoff) * S; //Units of A-1
						S = sqr(S);

						dedr = S * dedr + dSdd * epol;
						epol *= S;
					}
					dedr *= 2 * invdistij;
					dedx = dedr * dddx;
					dedy = dedr * dddy;
					dedz = dedr * dddz;

					epot_pol += epol;
					epot_pol += epol;
					epot_pol_cross += epol;
					epot_pol_cross += epol;

					GB_atom_param[i].dedi.x += dedx;
					GB_atom_param[i].dedi.y += dedy;
					GB_atom_param[i].dedi.z += dedz;
					GB_atom_param[j].dedi.x -= dedx;
					GB_atom_param[j].dedi.y -= dedy;
					GB_atom_param[j].dedi.z -= dedz;
				}

			}
			GB_atom_param[i].deda += dedai;

			epol = qi * wspace.atom[i].charge / ai;
			epot_pol += epol; // add up energies
			epot_pol_self += epol;

			GB_atom_param[i].deda += -epol / ai;
		}




		//add electrostatic energy to the total
		wspace.ene.epot += epot_pol;


		// ---------------------------------------------------------------------------
		// Now add the second part of the derivative
		int Type;
		double gpi;
		double invP5;
		double sqrtCCF, CCF, dCCF;
		double rdivR;
		double Rj, Ri;
		double Vj;
		double r6;
		double theta, costheta, sintheta;
		invP5 = 1 / P5;

		for(i = 0; i < wspace.atom.size(); i++) {
			Type = wspace.atom[i].FFType;
			Ri = GBtype[Type].radius;

			gpi = -sqr(GB_atom_param[i].bornradius) / -166.0;

			//all neighbors
			for(nj = 0; nj < fnbor[i].n; nj++) {
				j = NList32Bit_Index(fnbor[i].i[nj]);
				dddx = (wspace.cur.atom[j].p.x - wspace.cur.atom[i].p.x); // these are incomplete (need div by Dist_ij)
				dddy = (wspace.cur.atom[j].p.y - wspace.cur.atom[i].p.y);
				dddz = (wspace.cur.atom[j].p.z - wspace.cur.atom[i].p.z);
				srij = sqr(dddx) + sqr(dddy) + sqr(dddz);
				
				if(srij > sqr(GbsaStillCutoff))
					continue;
				//if( j > i ) continue;
				if(NList32Bit_BondOrder(fnbor[i].i[nj]) <= BondOrder_1_3_Pair)
					continue;// only get non-bonded and 1-4 neighbors
				
				Rj = GBtype[wspace.atom[j].FFType].radius;


				Vj = GB_atom_param_still[j].Vj;

				r6 = srij * srij * srij;
				rdivR = srij / sqr((Ri + Rj));

				if(rdivR > invP5) { // simple case, CCF = 1.0
					CCF = 1.0;
					dCCF = 0.0;
				} else {
					theta = rdivR * P5 * Maths::MathConst::PI;
					sintheta = sin(theta);
					costheta = cos(theta);
					sqrtCCF = 0.5 - 0.5 * costheta;
					CCF = sqr(sqrtCCF);
					dCCF = 2.0 * sqrtCCF * sintheta * theta;
				}


				dedr = gpi * GB_atom_param[i].deda * P4 * Vj * (4.0 * CCF - dCCF) / r6;

				dedx = dedr * dddx;
				dedy = dedr * dddy;
				dedz = dedr * dddz;

				GB_atom_param[i].dedi.x += dedx;
				GB_atom_param[i].dedi.y += dedy;
				GB_atom_param[i].dedi.z += dedz;
				GB_atom_param[j].dedi.x -= dedx;
				GB_atom_param[j].dedi.y -= dedy;
				GB_atom_param[j].dedi.z -= dedz;

			}
		}

		for(i = 0; i < atoms; i++) {
			force.setTo(GB_atom_param[i].dedi);
			force.mul(-1 / (PhysicsConst::Angstrom));
			wspace.cur.atom[i].f.add(force);
		}
		epot += epot_pol;
	}














	// -------------------------------------------------------------------------------
	// Class: FF_GeneralizedBorn_Still
	// Function: calcSelfBornEnergyDerivatives();
	// -------------------------------------------------------------------------------
	// Parameters: none
	//
	// this calculates the derivative of the self born terms wth respect to
	// the atom movement. It relies of the derivatives of the inverse born
	// radii as calculated by FF_GeneralizedBorn_Still::calcBornRadii_PairwiseApprox_PlusSelfDeriv();
	// it also assumes that the rest of the derivatives (the cross terms)
	// are calculated *before* !



	void FF_GeneralizedBorn_Still::calcForcesIncludingVacuo()
	{
		// Setup proxies to workspace to make code more readable.
		WorkSpace& wspace = getWSpace();
		size_t atoms = wspace.atom.size();              // number of atoms in workspace
		const ParticleStore& atomparam = wspace.atom; // atom parameter array
		SnapShotAtom *atom = wspace.cur.atom;           // atom coordinate array

		int i, j, nj;
		int typej;

		double dedr, dedalpha;
		double dedx, dedy, dedz;

		double dddx, dddy, dddz;

		const NeighbourData *fnbor = wspace.nlist().getData();

		double sqrcutoff = sqr(Cutoff);

		double Swidth = Cutoff - InnerCutoff;
		double invSwidth = 1.0 / Swidth;
		double S;
		double dSdd;

		double vdwSwidth = VdwCutoff - VdwInnerCutoff;
		double vdwinvSwidth = 1.0 / vdwSwidth;
		double vdwS;
		double vdwdSdd;

		double Dist_ij, invdistij;
		double qi, ai, aj, aiaj, srij;
		double uGB, fGB, lGB, mGB, tGB;
		double premul;
		double qiqj;
		double pair_elec14scaling = 1.0;

		double epol;
		double dedai;

		double atomix, atomiy, atomiz, force_magnitude;
		dvector fv, force;
		double vdw_force = 0, // individual force magnitudes
			elec_force = 0;
		double vdw_potential = 0, // individual potential energy contributions
			elec_potential = 0;

		double pair_vdw14scaling = 1.0;
		double radiusij, epsilon, A, B;
		double invdielectric = 1 / DielectricSolute;
		double invsolventdielectric = 1 / DielectricSolvent;

		double atomi_radius;
		double atomi_epsilon;
		dvector atomi_p;

		wspace.ene.epot_vdw = 0;
		wspace.ene.epot_elec = 0;

		epot = 0;
		epot_pol = 0.0;
		epot_pol_self = 0.0;
		epot_pol_cross = 0.0;

		atoms = wspace.atom.size();

		for(i = 0; i < atoms; i++)
			GB_atom_param[i].position.setTo(atom[i].p);

		premul = bornpremul;

		// set all arrays to 0
		for(i = 0; i < atoms; i++) 
		{
			GB_atom_param[i].dedi.setTo(0, 0, 0);
			GB_atom_param[i].deda = 0.0;
		}

		int pairs=0;

		for(i = 0; i < atoms; i++) 
		{
			ai = GB_atom_param[i].bornradius;
			atomix = atom[i].p.x;
			atomiy = atom[i].p.y;
			atomiz = atom[i].p.z;
			dedai = 0;

			atomi_radius = atomparam[i].radius;
			atomi_epsilon = atomparam[i].epsilon;
			qi = atomparam[i].charge;

			// do neighbors
			for(nj = 0; nj < fnbor[i].n; nj++) 
			{
				j = NList32Bit_Index(fnbor[i].i[nj]);
				if(j > i) break; // ignore shadow neighbors

				typej = NList32Bit_BondOrder(fnbor[i].i[nj]);
				
				aj = GB_atom_param[j].bornradius;
				qiqj = qi * GB_atom_param[j].charge;

				dddx = (double)(atomix - GB_atom_param[j].position.x); // these are incomplete (need div by Dist_ij)
				dddy = (double)(atomiy - GB_atom_param[j].position.y);
				dddz = (double)(atomiz - GB_atom_param[j].position.z);
				srij = sqr(dddx) + sqr(dddy) + sqr(dddz);
				if(srij>sqrcutoff) continue;
				Dist_ij = sqrt(srij);

				force_magnitude= 0; // total force/total pot. energy = 0

				invdistij = 1 / Dist_ij; // inverse distance - universally required

				aiaj = ai * aj;

				uGB = -srij * 0.25 / aiaj;

				pair_elec14scaling = 1.0;
				if(NList32Bit_BondOrder(fnbor[i].i[nj]) == BondOrder_1_4_Pair)
					pair_elec14scaling = Elec14Scaling;

				// approximate GB and vacuo electrostatic
				elec_potential = 0;
				elec_force = 0;
					
				if((uGB > ExpApproxThreshold)) 
				{
					tGB = exp(uGB);

					mGB = tGB * aiaj;
					lGB = 1/(srij + mGB);

					fGB = sqrt(lGB);
					epol = premul * qiqj * fGB;

					dedr = -epol * (Dist_ij - 0.25 * Dist_ij * tGB) * lGB;
					dedalpha = -epol * tGB * (1.0 - uGB) * lGB;
					/*
					lGB = srij + mGB;

					fGB = sqrt(lGB);
					epol = premul * qiqj / fGB;

					dedr = -epol * (Dist_ij - 0.25 * Dist_ij * tGB) / lGB;
					dedalpha = -epol * tGB * (1.0 - uGB) / lGB;

					*/
					if(typej >= 3) 
					{ 
						// only do vdw/elec for non-bonded and 1-4 neighbors
						elec_potential = (double) PhysicsConst::econv_joule *invdielectric * pair_elec14scaling * qiqj * invdistij;
						elec_force = -elec_potential * invdistij;
						pairs++;
					}
					// Apply switching function (S)
					if(Dist_ij >= InnerCutoff) 
					{
						S = (1.0 - sqr(invSwidth * (Dist_ij - InnerCutoff))); //Dimension less
						dSdd = -4.0 * sqr(invSwidth) * (Dist_ij - InnerCutoff) * S; //Units of A-1
						S = sqr(S);

						elec_force = S * elec_force + elec_potential * dSdd;
						elec_potential *= S;
						dedr = S * dedr + dSdd * epol;
						dedalpha *= S;
						epol *= S;
					}

					dedr *= 2 * invdistij;
					dedx = dedr * dddx;
					dedy = dedr * dddy;
					dedz = dedr * dddz;

					epot_pol += epol;
					epot_pol += epol;
					epot_pol_cross += epol;
					epot_pol_cross += epol;
					wspace.atom[i].epot += epol;
					wspace.atom[j].epot += epol;

					GB_atom_param[i].dedi.x += dedx;
					GB_atom_param[i].dedi.y += dedy;
					GB_atom_param[i].dedi.z += dedz;
					GB_atom_param[j].dedi.x -= dedx;
					GB_atom_param[j].dedi.y -= dedy;
					GB_atom_param[j].dedi.z -= dedz;

					dedai += dedalpha * aj;
					GB_atom_param[j].deda += dedalpha * ai;

					
				} 
				else 
				{
					// water Dielectric electrostatics (long distance)
					if(typej >= 3) 
					{ 
						// only do vdw/elec for non-bonded and 1-4 neighbors
						elec_potential = PhysicsConst::econv_joule * invsolventdielectric * pair_elec14scaling * qiqj * invdistij;
						elec_force = -elec_potential * invdistij;

						// Apply switching function (S)
						if(Dist_ij >= InnerCutoff) 
						{
							S = (1.0 - sqr(invSwidth * (Dist_ij - InnerCutoff))); //Dimension less
							dSdd = -4.0 * sqr(invSwidth) * (Dist_ij - InnerCutoff) * S; //Units of A-1
							S = sqr(S);
							elec_force = S * elec_force + elec_potential * dSdd;
							elec_potential *= S;
						}
					}
				}

				pair_vdw14scaling = 1.0;
				if(typej >= 3) { // only do vdw/elec for non-bonded and 1-4 neighbors
					if(typej == 3)
						pair_vdw14scaling = Vdw14Scaling;

					vdw_force = 0;
					if(Dist_ij < VdwCutoff) {
						// Van der Waals Force
						radiusij = sqr(invdistij * (atomi_radius + GB_atom_param[j].radiusij));
						epsilon = pair_vdw14scaling * atomi_epsilon * GB_atom_param[j].epsilon;

						//invd6 = cube()*cube(invdistij); // 1/d^6 needed for vdw
						B = cube(radiusij);
						A = sqr(B);

						vdw_potential = (2.0 * epsilon) * (0.5 * A - B);
						vdw_force = -(12.0 * epsilon * invdistij) * (A - B);


						if(Dist_ij > VdwInnerCutoff) {
							vdwS = (1.0 - sqr(vdwinvSwidth * (Dist_ij - VdwInnerCutoff))); // Dimension less
							vdwdSdd = -4.0 * sqr(vdwinvSwidth) * (Dist_ij - VdwInnerCutoff) * vdwS; // Units of per meter (not per angstrom!)
							vdwS = sqr(vdwS);
							vdw_force = vdwS * vdw_force + vdw_potential * vdwdSdd;

							vdw_potential *= vdwS;
						}
						// printf(" %lf %lf \n",VdwCutoff, VdwInnerCutoff);

						// add up potentials -----------------------------

						wspace.ene.epot += vdw_potential;
						wspace.ene.epot_vdw += vdw_potential;
						wspace.atom[i].epot += vdw_potential * 0.5;
						wspace.atom[j].epot += vdw_potential * 0.5;
					}
					wspace.ene.epot += elec_potential;
					wspace.ene.epot_elec += elec_potential;
					wspace.atom[i].epot += elec_potential * 0.5;
					wspace.atom[j].epot += elec_potential * 0.5;

					// now apply the force
					force_magnitude= 0;
					force_magnitude+= vdw_force;
					force_magnitude+= elec_force;
					force_magnitude*= invdistij;
					dedx = force_magnitude* dddx;
					dedy = force_magnitude* dddy;
					dedz = force_magnitude* dddz;

					GB_atom_param[i].dedi.x += dedx;
					GB_atom_param[i].dedi.y += dedy;
					GB_atom_param[i].dedi.z += dedz;
					GB_atom_param[j].dedi.x -= dedx;
					GB_atom_param[j].dedi.y -= dedy;
					GB_atom_param[j].dedi.z -= dedz;
				}
			}
			GB_atom_param[i].deda += dedai;

			epol = premul * sqr(qi) / ai;
			epot_pol += epol; // add up energies
			epot_pol_self += epol;
			wspace.atom[i].epot += epol;

			GB_atom_param[i].deda += -epol / ai;
		}

		//add electrostatic energy to the total
		wspace.ene.epot += epot_pol;
		//add electrostatic energy to the energy components
		wspace.ene.epot_pol += epot_pol;
		wspace.ene.epot_pol_cross += epot_pol_cross;
		wspace.ene.epot_pol_self += epot_pol_self;
		//printf("pairs: %d \n",pairs);
		// ---------------------------------------------------------------------------
		// Now add the second part of the derivative
		int Type;
		double gpi;
		double invP5;
		double sqrtCCF, CCF, dCCF;
		double rdivR;
		double Rj, Ri;
		double Vj;
		double r6;
		double theta, costheta, sintheta;
		invP5 = 1 / P5;

		for(i = 0; i < wspace.atom.size(); i++) {
			Type = wspace.atom[i].FFType;
			Ri = GBtype[Type].radius;

			gpi = -sqr(GB_atom_param[i].bornradius) / -166.0;

			//all neighbors
			for(nj = 0; nj < fnbor[i].n; nj++) {
				j = NList32Bit_Index(fnbor[i].i[nj]);
				dddx = (wspace.cur.atom[j].p.x - wspace.cur.atom[i].p.x); // these are incomplete (need div by Dist_ij)
				dddy = (wspace.cur.atom[j].p.y - wspace.cur.atom[i].p.y);
				dddz = (wspace.cur.atom[j].p.z - wspace.cur.atom[i].p.z);
				srij = sqr(dddx) + sqr(dddy) + sqr(dddz);
				
				if(srij > sqr(GbsaStillCutoff))
					continue;
				//if( j >= i) continue;
				if(NList32Bit_BondOrder(fnbor[i].i[nj]) <= BondOrder_1_3_Pair)
					continue;// only get non-bonded and 1-4 neighbors
				
				// but treat shadow and normal neighbors the same
				Rj = GBtype[wspace.atom[j].FFType].radius;


				Vj = GB_atom_param_still[j].Vj;

				//srij = sqr(rij);
				r6 = srij * srij * srij;
				rdivR = srij / sqr((Ri + Rj));

				if(rdivR > invP5) { // simple case, CCF = 1.0
					CCF = 1.0;
					dCCF = 0.0;
				} else {
					theta = rdivR * P5 * Maths::MathConst::PI;
					sintheta = sin(theta);
					costheta = cos(theta);
					sqrtCCF = 0.5 - 0.5 * costheta;
					CCF = sqr(sqrtCCF);
					dCCF = 2.0 * sqrtCCF * sintheta * theta;
				}


				dedr = gpi * GB_atom_param[i].deda * P4 * Vj * (4.0 * CCF - dCCF) / r6;

				dedx = dedr * dddx;
				dedy = dedr * dddy;
				dedz = dedr * dddz;

				GB_atom_param[i].dedi.x += dedx;
				GB_atom_param[i].dedi.y += dedy;
				GB_atom_param[i].dedi.z += dedz;
				GB_atom_param[j].dedi.x -= dedx;
				GB_atom_param[j].dedi.y -= dedy;
				GB_atom_param[j].dedi.z -= dedz;

			}
		}

		for(i = 0; i < atoms; i++) {
			force.setTo(GB_atom_param[i].dedi);
			force.mul((double) - 1.0 / (double) (PhysicsConst::Angstrom));
			wspace.cur.atom[i].f.add(force);
		}
    epot_elec = wspace.ene.epot_elec;
    epot_vdw  = wspace.ene.epot_vdw; 	

	}
} // namespace Physics

