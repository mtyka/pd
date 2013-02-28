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

#include "forcefields/ffparam.h"
#include "workspace/workspace.h"
#include "workspace/neighbourlist.h"

// self header include
#include "lcpo.h"

// namespace includes
using namespace Maths;

namespace Physics
{
	void FF_SASA_LCPO::setup()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		wspace.nlist().requestCutoff(7.0);

		// Alloc
		SASAatom.clear();
		SASAtype.clear();
		Nlist.clear();
		SASAatom.resize(wspace.atom.size());
		//SASAtype.resize(wspace.ffps().AtomType.size() * 4);
		Nlist.resize(wspace.atom.size());

		int attype;
		int neighbors;

		// read algorithmic parameters
		// this one must come first, since it 'creates' the SASA types, all other ones have to be later and
		// relate to the present SASA types !
		if(readLCPOSection() != 0)  throw(ProcedureException("Error reading Solvation parameters")); 

		if(!cmpstring(ASPsection_name,""))
		{
			if(readASPSection(ASPsection_name) != 0)  throw(ProcedureException("Error reading Solvation parameters"));
		} 
		else 
		{
			for(size_t i = 0; i < SASAtype.size(); i++)
			{
				SASAtype[i].sigma = GlobalASP / PhysicsConst::Na / PhysicsConst::J2kcal;
			}
		}

		for(size_t i = 0; i < wspace.atom.size(); i++) 
		{
			attype = wspace.atom[i].FFType; // get attype
			neighbors = 0; // count neighbors
			for(size_t j = 0; j < wspace.atom[i].cov12atom.size(); j++) 
			{
				if(wspace.atom[wspace.atom[i].cov12atom[j].i].Z != 1)
					neighbors++;
			}

			size_t j; // required outside the loop
			for(j = 0; j < SASAtype.size(); j++) 
			{
				if((SASAtype[j].attype == attype)) 
				{
					if((SASAtype[j].neighbors < 0) || (SASAtype[j].neighbors == neighbors))
					{
						break;
					}
				}
			}

			if(j >= SASAtype.size()) 
			{
				printf("ERROR: No appropriate SASAtype for Atom %d found (%s %d %s )\n", i,
					wspace.atom[i].parentl3name.c_str(),
					wspace.atom[i].ir,
					wspace.atom[i].pdbname.c_str());
				 throw(ProcedureException("Error reading Solvation parameters"));	
			}

			// EUGH -> memcpy((void *) &SASAatom[i], (void *) &SASAtype[j], sizeof(SASA_Atom));
			SASAatom[i] = SASAtype[j];
			SASAatom[i].radius += 1.4;
			// ensure we ignore hydrogens (Z==1);
			if(wspace.atom[i].Z == 1)
				SASAatom[i].use = 0;
			else
				SASAatom[i].use = 1;
		}


		/*
		for(i=0;i<wspace.atom.size();i++){
		printf(" %4d %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %2d %2d \n",
		SASAatom[i].attype, // link to the atomt ype
		SASAatom[i].P1,
		SASAatom[i].P2,
		SASAatom[i].P3,
		SASAatom[i].P4, // LCPO Parameters
		SASAatom[i].radius, // Vdw Radius + 1.4
		SASAatom[i].sigma, // atomic solvation parameter
		SASAatom[i].neighbors, // Nr of neighbors |???
		SASAatom[i].use); // include in calculation
		}

		*/

		needsetup = false;
		Active = true;
	}

	void FF_SASA_LCPO::info() const
	{ 
		// prints a little block of parameter information
		printf(" --SASA Forcefield ------------------ \n");
		if(!Active)
			printf(" -----------------------<INACTIVE>---\n");
		if(cmpstring(ASPsection_name,"")){
			printf(" Sigma: %8.6lf kcal/mol\n", GlobalASP );
		}else{
			printf(" FF Section: %s\n",ASPsection_name.c_str() );
		}
		printf("\n");
	}

	void FF_SASA_LCPO::infoLine() const 
	{ 
		// prints a line of current energies
		printf("\t%6.1lf", this->epot_cav * PhysicsConst::J2kcal * PhysicsConst::Na);
	}

	void FF_SASA_LCPO::infoLineHeader() const 
	{ 
		// prints the headers for the above function
		printf("\t%6s", "Esurf");
	}

	// This displays a table of energies in a pairwise fashion - not for use
	// use in MD but for the user's single energy evaluaions
	void FF_SASA_LCPO::calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level)
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		if(level == Summary) {
			calcEnergies();
			printf(" SASA Energy:   %10.3lf kcal/mol (%.1lf A^2)\n",
				epot_cav * PhysicsConst::J2kcal * PhysicsConst::Na,
				totalSASA );
			return;
		}

		const NeighbourData *fnbor = wspace.nlist().getData();
		calcLCPOSasaForces(true);
		printf(" Nr Neigbrs SASArad SASA ASP(cal) Ecav alpha Vj Term123 Charge self offdiag epol \n\n");
		for(int i = 0; i < wspace.atom.size(); i++) {
			printf("SoluteAtom: %4d %6d %6.3lf %6.3lf %5.2lf %6.4lf\n",
				i,
				fnbor[i].n,
				SASAatom[i].radius,
				SASAatom[i].SASA,
				SASAatom[i].sigma * PhysicsConst::Na * PhysicsConst::J2kcal * 1000,
				SASAatom[i].SASA * SASAatom[i].sigma * PhysicsConst::Na * PhysicsConst::J2kcal);
		}
		wspace.ene.epot_surf += epot_cav;
		epot = epot_cav;
	}

	void FF_SASA_LCPO::calcEnergies()
	{
		calcLCPOSasaEnergies(); // calculates SASA & adds energies to wspace
		getWSpace().ene.epot_surf += epot_cav;
		epot = epot_cav;
	}

	void FF_SASA_LCPO::calcForces()
	{
		WorkSpace& wspace = getWSpace();
		if((wspace.Step % UpdateSasa) == 0)
		{
			calcLCPOSasaForces(true); // calculates SASA & adds energies to wspace
		}
		else
		{
			calcLCPOSasaForces(false); // adds previous energies to wspace
		}
		wspace.ene.epot_surf += epot_cav;
		epot = epot_cav;
	}

	void FF_SASA_LCPO::testDerivatives()
	{
		printf("----\n");
		calcLCPOSasaEnergies(); // calculates SASA & adds energies to wspace
		printf("----\n");
		calcLCPOSasaForces(true);
		printf("----\n");
		calcLCPOSasaForces_num();
		printf("----\n");
	}






	int FF_SASA_LCPO::readASPSection(const std::string &sectionname)
	{
		int done;
		double solvationParameter;

		int j,line, attype;
		int errorstatus = 0;

		Section section;
		if(getWSpace().ffps().getSection(sectionname,section)!=0){
			printf("ERROR: Cannot find section %s in forcefield parameters\n",
				sectionname.c_str());
			return -1;
		}

		//Start reading parameter file

		line = section.lineoffset;
		for(unsigned i=0;i<section.sectionline.size();i++){
			std::string linestring = section.sectionline[i];
			removecomments(linestring,"#\12\15");
			line++;
			std::vector<std::string> token;
			token = chopstr(linestring," \12\15\t");
			if(token.size() <= 0) continue; // empty line
			std::string command = token[0];

			if(cmpstring(command, "SOLVATIONTYPE")) {

				if(token.size() < 3) {
					printf
						("SYNTAX ERROR (line %d): Insufficient parameters for SOLVATIONTYPE (either 2 or 7 parameters required)\n",
						line);
					errorstatus = 1;
					continue;
				}
				attype = getWSpace().ffps().findAtomType(token[1]); // find attype number
				if(attype < 0) { // in case we didn't find it
					printf("SYNTAX ERROR (line %d): Atom name %s unknown \n", line, token[1].c_str());
					errorstatus = 1;
					continue;
				}

				if(str2double(token[2],solvationParameter)!=0){ // assign parameters
					printf("ERROR: Numerical value expected for parameter 2 after SOLVATIONTYPE\n");
					errorstatus = 1;
					continue;
				}

				done = 0;
				for(j = 0; j < SASAtype.size(); j++) 
				{
					if(SASAtype[j].attype == attype) 
					{ 
						// found it in the curretn sasa types
						SASAtype[j].sigma = solvationParameter * PhysicsConst::kcal2J * PhysicsConst::Na;
						done = 1;
					}
				}
				if(done == 0) { // in case we didn't find it
					printf("SYNTAX ERROR (line %d): Atom name %s unknown \n", line, token[1].c_str());
					errorstatus = 1;
					continue;
				}


			} else {
				printf("ERROR: Identifier %s not known in %s section \n",
					command.c_str(), sectionname.c_str());
				errorstatus = 1;
			}
		}


		// Error handling

		if(errorstatus == 0){
			printf("Finished reading ASP solvation parameters - no errors\n");
		}else{
			printf("Errors during reading of Eisenberg et al. solvation parameters - check script \n");
		}

		return errorstatus;
	}


	int FF_SASA_LCPO::readLCPOSection()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		double P1, P2, P3, P4, radius;
		int line,neighbors,attype;
		int errorstatus = 0;
		std::string sectionname = "LCPOSASAPARAMETERS";

		Section section;
		if(wspace.ffps().getSection(sectionname,section)!=0){
			printf("ERROR: Cannot find section %s in forcefield parameters\n",
				sectionname.c_str());
			return -1;
		}

		//Start reading parameter file

		line = section.lineoffset;
		for(unsigned i=0;i<section.sectionline.size();i++)
		{
			std::string linestring = section.sectionline[i];
			removecomments(linestring,"#\12\15");
			line++;
			std::vector<std::string> token;
			token = chopstr(linestring," \12\15\t");
			if(token.size() <= 0) continue; // empty line
			std::string command = token[0];

			if(cmpstring(command, "LCPOSASATYPE")) { // TYPE syntax reading/checking of SOLVATIONTYPE

				if(token.size() < 8) {
					printf("SYNTAX ERROR (line %d): Insufficient parameters for LCPOSASATYPE \n", line);
					errorstatus = 1;
					continue;
				}

				attype = wspace.ffps().findAtomType(token[1]); // find attype number
				if(attype < 0) { // in case we didn't find it
					printf("SYNTAX ERROR (line %d): Atom name %s unknown \n", line, token[1].c_str());
					errorstatus = 1;
					continue;
				}

				if(str2int(token[2],neighbors)!=0){ // assign parameters
					printf("ERROR: Numerical value expected for number of neighbors after LCPOSASATYPE\n");
					errorstatus = 1;
					continue;
				}
				if(str2double(token[3],radius)!=0){ // assign parameters
					printf("ERROR: Numerical value expected for parameter radius after LCPOSASATYPE\n");
					errorstatus = 1;
					continue;
				}
				if(str2double(token[4],P1)!=0){ // assign parameters
					printf("ERROR: Numerical value expected for parameter P1 after LCPOSASATYPE\n");
					errorstatus = 1;
					continue;
				}
				if(str2double(token[5],P2)!=0){ // assign parameters
					printf("ERROR: Numerical value expected for parameter P2 after LCPOSASATYPE\n");
					errorstatus = 1;
					continue;
				}
				if(str2double(token[6],P3)!=0){ // assign parameters
					printf("ERROR: Numerical value expected for parameter P3 after LCPOSASATYPE\n");
					errorstatus = 1;
					continue;
				}
				if(str2double(token[7],P4)!=0){ // assign parameters
					printf("ERROR: Numerical value expected for parameter P4 after LCPOSASATYPE\n");
					errorstatus = 1;
					continue;
				}

				SASA_Atom a;
				a.radius = radius;
				a.attype = attype;
				a.neighbors = neighbors;
				a.P1 = P1;
				a.P2 = P2;
				a.P3 = P3;
				a.P4 = P4;
				SASAtype.push_back(a);

			} 
			else 
			{
				printf("ERROR: Identifier %s not known in %s section \n",
					command.c_str(), sectionname.c_str());
				errorstatus = 1;
			}
		}


		if(errorstatus == 0){
			printf("Finished reading LCPO (SASA) parameter file - no errors\n");
		}else{
			printf("Errors during reading of LCPO (SASA) parameter file - check script \n");
		}

		/*
		for(i=0;i<nSASAtypes;i++){
		if(SASAtype[i].radius < 0) printf("SASAtype %d -- invalid\n",i); // invalid Type
		else printf("SASAtype %d %d %d %lf %lf %lf %lf %lf \n",
		i,
		SASAtype[i].attype,
		SASAtype[i].neighbors,
		SASAtype[i].radius,
		SASAtype[i].P1,
		SASAtype[i].P2,
		SASAtype[i].P3,
		SASAtype[i].P4);

		}
		*/

		return errorstatus;
	}

	int FF_SASA_LCPO::calcLCPOSasaEnergies(bool recalc_nlist)
	{
		// a fast inline implementation of lcpoSASA to work with
		// existing neighbor lists etc

		int i, j, k, ia, ja, ka;
		double SAi, Si, Si2, Si3, Si4, Si2_t, Si3_t;
		int count = 0;

		// Proxies
		WorkSpace& wspace = getWSpace();
		int atoms = wspace.atom.size();
		const NeighbourData *fnbor = wspace.nlist().getData();
		double Dist_ij;

		totalSASA = 0;
		epot_cav = 0;
int ncount=0;
		if(recalc_nlist){
			// create neighbour lists
			count = 0;
			for(i = 0; i < atoms; i++)
				Nlist[i].size = 0;

			for(i = 0; i < atoms; i++) {
				if(SASAatom[i].use == 0) {
					count += atoms;
					continue;
				}
				count += i + 1;

				//run over all neighbors
				for(int nj = 0; nj < fnbor[i].n; nj++) {
					j = NList32Bit_Index(fnbor[i].i[nj]);
					if(j >= i) break;
					//if(NList32Bit_BondOrder(fnbor[i].i[nj]) <= 3) continue;
					// old:
					//if(fnbor[i].Type[nj] >= 128)
					//	continue; // only get direct neighbors, not shadow neighbors
					//j = fnbor[i].i[nj];
					
					ncount++;
					if(SASAatom[j].use == 0) {
						count++;
						continue;
					}
					Dist_ij = dist(wspace.cur.atom[i].p,wspace.cur.atom[j].p);

					// for(j = i+1;j<atoms;j++){
					// if(SASAatom[j].use==0){ count++; continue; }
					// Dist_ij = dmatrix[count];

					if(Dist_ij < (SASAatom[i].radius + SASAatom[j].radius)) {
						if(Nlist[i].size<Nlistmax){
							Nlist[i].i[Nlist[i].size] = j;
							Nlist[i].s[Nlist[i].size] = sphereSurfaceOverlap(SASAatom[i].radius, SASAatom[j].radius, Dist_ij);
							if(Nlist[i].s[Nlist[i].size] < 0)
								Nlist[i].s[Nlist[i].size] = 0;

							Nlist[i].size++;
						}
						if(Nlist[j].size<Nlistmax){
							Nlist[j].i[Nlist[j].size] = i;
							Nlist[j].s[Nlist[j].size] = sphereSurfaceOverlap(SASAatom[j].radius, SASAatom[i].radius, Dist_ij);
							if(Nlist[j].s[Nlist[j].size] < 0)
								Nlist[j].s[Nlist[j].size] = 0;
							Nlist[j].size++;
						}
					}
					count++;
				}
			}

		}else{
			for(i = 0; i < atoms; i++) {
				for(j = 0; j < Nlist[i].size; j++) {
					ja = Nlist[i].i[j]; // get atom number
					Dist_ij = wspace.cur.atom[i].p.dist(wspace.cur.atom[ja].p);
					Nlist[i].s[j] = 
						sphereSurfaceOverlap(SASAatom[i].radius, 
						SASAatom[ja].radius, Dist_ij);

				}
			}

		}

		for(i = 0; i < atoms; i++) {
			if(SASAatom[i].use == 0) {
				SASAatom[i].SASA = 0;
				continue; // ignore radius 0 atoms
			}

			ia = i;
			Si = sphereSurfaceArea(SASAatom[ia].radius);

			Si2 = 0;
			Si3 = 0;
			Si4 = 0;

			for(j = 0; j < Nlist[i].size; j++) {
				ja = Nlist[ia].i[j]; // get atom number

				Si2_t = Nlist[ia].s[j];
				Si3_t = 0;
				for(k = 0; k < Nlist[ja].size; k++) {
					ka = Nlist[ja].i[k];
					if(ka == ia)
						continue; // exclude if ia and ka are the same atom

					//check if those two (ka and ia) actually overlap
					if((SASAatom[ia].radius + SASAatom[ka].radius) < 
						wspace.cur.atom[ia].p.dist(wspace.cur.atom[ka].p)
						//dmatrix[sqrmat(ia,ka,atoms)]
						)
						continue;

					Si3_t += Nlist[ja].s[k];
				}

				Si4 += (Si3_t * Si2_t);
				Si3 += Si3_t;
				Si2 += Si2_t;
			}

			SAi = SASAatom[i].P1 * Si +
				SASAatom[i].P2 * Si2 +
				SASAatom[i].P3 * Si3 +
				SASAatom[i].P4 * Si4;

			if(SAi < 0)
				SAi = 0;
			// printf("SASA = %8.2lf -- Si %6.2lf Si2 %6.2lf Si3 %6.2lf Si4 %6.2lf \n",SAi,Si,Si2,Si3,Si4);

			SASAatom[i].SASA = SAi;

			// add things to wspace
			epot_cav += (double) SASAatom[i].SASA * SASAatom[i].sigma;

			totalSASA += SAi;
		}

		wspace.ene.epot += epot_cav;

		return 0;
	}







	int FF_SASA_LCPO::calcLCPOSasaForces_num(){// a slow version which calculates the derivatives of the SASA of each
		int k;
		double dc = 0.001; // displacement in Angstrom
		dvector atompos;

		// Proxies
		WorkSpace& wspace = getWSpace();
		int atoms = wspace.atom.size();
		dvector *depotdc = new dvector[atoms]; // store the changed epot for each atom move
		double epot_0;

		// Ok, now calculate each atom's SASA when nothing has moved

		calcLCPOSasaEnergies(true); 
		epot_0 = epot_cav;
		printf("\n");
		for(k = 0; k < atoms; k++) {
			depotdc[k].setTo(0,0,0);

			if(SASAatom[k].use == 0) continue;

			atompos.setTo(wspace.cur.atom[k].p);  // save the original position

			wspace.cur.atom[k].p.x += dc; // add an x displacement
			wspace.nlist().calcNewList();
			calcLCPOSasaEnergies(false); 
			depotdc[k].x = (epot_0 - epot_cav)/(dc*PhysicsConst::Angstrom);
			wspace.cur.atom[k].p.setTo(atompos); // reset the moved atom

			wspace.cur.atom[k].p.y += dc;
			wspace.nlist().calcNewList();
			calcLCPOSasaEnergies(false); 
			depotdc[k].y = (epot_0 - epot_cav)/(dc*PhysicsConst::Angstrom);
			wspace.cur.atom[k].p.setTo(atompos); // reset the moved atom

			wspace.cur.atom[k].p.z += dc;
			wspace.nlist().calcNewList();
			calcLCPOSasaEnergies(false); 
			depotdc[k].z = (epot_0 - epot_cav)/(dc*PhysicsConst::Angstrom);
			wspace.cur.atom[k].p.setTo(atompos); // reset the moved atom

			if(SASAatom[k].use == 0) continue;
			printf("numerical: %d %e %e %e \n",k, depotdc[k].x,depotdc[k].y,depotdc[k].z );
			wspace.cur.atom[k].f.add(depotdc[k]);
		}

		calcLCPOSasaEnergies(false); 

		delete[]depotdc;
		return 0;
	}





	int FF_SASA_LCPO::calcLCPOSasaForces(bool dofullcalc)
	{
		// a fast inline implementation of lcpoSASA to work with
		// existing neighbor lists etc
		// this version also calculates the derivatives of the SASA of each
		// atom with respect to that atoms' position
		int i, j, k, ia, ja, ka;
		double SAi, Si, Si2, Si3, Si4;
		double Si2_t;
		double Si3_t;
		double Si3ik_t;
		dvector forcem;

		// Proxies
		WorkSpace& wspace = getWSpace();
		int atoms = wspace.atom.size();
		const NeighbourData *fnbor = wspace.nlist().getData();

		SAi = 0;
		if(!dofullcalc) 
		{
			// unless instructured to do a full calculation
			// just use forces & energies from last Step
			epot_cav = 0.0;
			for(i = 0; i < atoms; i++) 
			{
				epot_cav += (double) SASAatom[i].SASA * SASAatom[i].sigma;
				wspace.atom[i].epot += (double) SASAatom[i].SASA * SASAatom[i].sigma;
				forcem.setTo(SASAatom[i].SASAderiv);
				forcem.mul(1 / (PhysicsConst::Angstrom));
				wspace.cur.atom[i].f.add(forcem);

				// sum up total SASA
				totalSASA += SAi;
			}

			wspace.ene.epot += epot_cav;
			return 0;
		}

		double	sij, sji, sjk;
		int			count = 0;
		double	Dist_ij, invdistij, distik, invdistik;
		double	ri, rj, rk;

		dvector Atom_i, Atom_j, atomk;

		dvector dAijdc_2t, dAijdc_2;
		dvector dAijdc_4t, dAijdc_4;

		dvector dSASA_2_neigh_dc;
		dvector dSASA_3_neigh_dc;
		dvector dSASA_4_neigh_dc;
		dvector dSASA_3_neigh_dc2;
		dvector dSASA_4_neigh_dc2;

		double	dAdd;
		double	dddx, dddy, dddz;
		double	ri2rj2divd2;

		member *Nlisti, *Nlistj;
		int Nlistsizei, Nlistsizej;

		// set the cavity/vdw potential to 0
		epot_cav = 0;
		totalSASA = 0;

		// create close contact neighbour lists
		count = 0;
		for(i = 0; i < atoms; i++)
			Nlist[i].size = 0;

		for(i = 0; i < atoms; i++) {
			if(SASAatom[i].use == 0) {
				count += atoms;
				continue;
			}
			count += i + 1;
			// load atom coords
			Atom_i.x = wspace.cur.atom[i].p.x;
			Atom_i.y = wspace.cur.atom[i].p.y;
			Atom_i.z = wspace.cur.atom[i].p.z;

			ri = SASAatom[i].radius;

			//run over all neighbors
			for(int nj = 0; nj < fnbor[i].n; nj++) {
				j = NList32Bit_Index(fnbor[i].i[nj]);
				if(j >= i) break;
				//if(fnbor[i].Type[nj] > 126)
				//	continue; // only get direct neighbors, not shadow neighbors
				//j = fnbor[i].i[nj];

				if(SASAatom[j].use == 0) {
					count++;
					continue;
				}
				Dist_ij = sqrdist(wspace.cur.atom[i].p,wspace.cur.atom[j].p);

				rj = SASAatom[j].radius;
				if(Dist_ij < sqr(ri + rj)) {
					Dist_ij = sqrt(Dist_ij);
					// load atom coords
					Atom_j.x = wspace.cur.atom[j].p.x;
					Atom_j.y = wspace.cur.atom[j].p.y;
					Atom_j.z = wspace.cur.atom[j].p.z;

					invdistij = 1 / Dist_ij;

					Nlisti = &Nlist[i];
					Nlistsizei = Nlisti->size;

					Nlistj = &Nlist[j];
					Nlistsizej = Nlistj->size;

					ri2rj2divd2 = 0.5 * (sqr(ri) - sqr(rj)) * invdistij;
					sij = Maths::MathConst::TwoPI * ri * (ri - Dist_ij * 0.5 - ri2rj2divd2);
					sji = Maths::MathConst::TwoPI * rj * (rj - Dist_ij * 0.5 + ri2rj2divd2);

					if(Nlist[i].size<Nlistmax){

						Nlisti->i[Nlistsizei] = j; // Load atom index
						Nlisti->s[Nlistsizei] = sij; // Load atom atom sphere overlap
						Nlisti->size++;
					}
					// -------------------------------

					if(Nlist[j].size<Nlistmax){
						Nlistj->i[Nlistsizej] = i; // Load atom index
						Nlistj->s[Nlistsizej] = sji; // Load atom atom sphere overlap
						Nlistj->size++;
					}
				}
			}
		}


		for(i = 0; i < atoms; i++) {

			Atom_i.setTo(wspace.cur.atom[i].p);

			//set derivative of this atom to 0
			SASAatom[i].SASAderiv.x = 0;
			SASAatom[i].SASAderiv.y = 0;
			SASAatom[i].SASAderiv.z = 0;

			if(SASAatom[i].use == 0) {
				SASAatom[i].SASA = 0;
				continue; // ignore radius 0 atoms
			}

			ia = i;
			ri = SASAatom[ia].radius;
			Si = sphereSurfaceArea(ri);

			Si2 = 0;
			Si3 = 0;
			Si4 = 0;

			dAijdc_2.setTo(0, 0, 0);
			dAijdc_4.setTo(0, 0, 0);

			for(j = 0; j < Nlist[i].size; j++) {
				ja = Nlist[ia].i[j]; // get atom number
				Atom_j.setTo(wspace.cur.atom[ja].p);
				Si2_t = Nlist[ia].s[j]; // get overlap with that neighbor
				Si3_t = 0;
				Si3ik_t = 0;
				Dist_ij = Atom_i.dist(Atom_j);

				rj = SASAatom[ja].radius;

				invdistij = 1 / Dist_ij;

				ri2rj2divd2 = 0.5 * (sqr(ri) - sqr(rj)) * invdistij;
				sji = Maths::MathConst::TwoPI * rj * (rj - Dist_ij * 0.5 + ri2rj2divd2);

				ri2rj2divd2 *= 2.0 * invdistij;

				dSASA_3_neigh_dc.setTo(0, 0, 0); // set third term of derivative to 0
				dSASA_4_neigh_dc.setTo(0, 0, 0); // set fourth of derivative to 0
				dSASA_3_neigh_dc2.setTo(0, 0, 0); // set third term of derivative to 0
				dSASA_4_neigh_dc2.setTo(0, 0, 0); // set fourth of derivative to 0

				//direct derivatives on i's neighbor j
				dAdd = Maths::MathConst::PI * rj * (-ri2rj2divd2 - 1);

				dddx = -(Atom_i.x - Atom_j.x) / Dist_ij;
				dddy = -(Atom_i.y - Atom_j.y) / Dist_ij;
				dddz = -(Atom_i.z - Atom_j.z) / Dist_ij;

				dSASA_2_neigh_dc.x = dAdd * dddx; // multiply with ddij / di
				dSASA_2_neigh_dc.y = dAdd * dddy;
				dSASA_2_neigh_dc.z = dAdd * dddz;

				// derivatives on i itself
				dAdd = Maths::MathConst::PI * ri * (ri2rj2divd2 - 1);
				dAijdc_2t.x = dAdd * dddx;
				dAijdc_2t.y = dAdd * dddy;
				dAijdc_2t.z = dAdd * dddz;

				for(k = 0; k < Nlist[ja].size; k++) {
					ka = Nlist[ja].i[k];

					atomk.setTo(wspace.cur.atom[ka].p);
					rk = SASAatom[ka].radius;

					// exclude if ia and ka are the same atom
					if(ka == ia)continue;

					// check if those two (ka and ja) actually overlap
					// i.e. k must be both neighbor of j as well as neighbor of i
					distik = Atom_i.sqrdist(atomk);

					if(sqr(ri + rk) < distik)	continue;
					distik = sqrt(distik);
					// indirect derivatives on neighbor j through k (j's neighbor)
					invdistik = 1 / distik;
					ri2rj2divd2 = (sqr(ri) - sqr(rk)) * sqr(invdistik);
					// calculate the derivative of the pairwise ik overlap with respect
					// to the ik distance
					dAdd = Maths::MathConst::PI * ri * (ri2rj2divd2 - 1);

					//calculate the derivative of distik with respect to atom position of ia
					dddx = -(Atom_i.x - atomk.x) * invdistik;
					dddy = -(Atom_i.y - atomk.y) * invdistik;
					dddz = -(Atom_i.z - atomk.z) * invdistik;

					// add up components of fourth term of derivative
					sjk = Nlist[ja].s[k];
					// add up components of third term of derivative
					dSASA_3_neigh_dc.x += dAdd * dddx; // multiply with ddij / di
					dSASA_3_neigh_dc.y += dAdd * dddy;
					dSASA_3_neigh_dc.z += dAdd * dddz;

					dAdd = Maths::MathConst::PI * rk * (-ri2rj2divd2 - 1);
					Si3ik_t += Maths::MathConst::PI * ri * (2*ri - distik  - distik*ri2rj2divd2);	

					// add up components of third term of derivative
					dSASA_3_neigh_dc2.x += dAdd * dddx; // multiply with ddij / di
					dSASA_3_neigh_dc2.y += dAdd * dddy;
					dSASA_3_neigh_dc2.z += dAdd * dddz;

					dSASA_4_neigh_dc2.x += sjk * dAdd * dddx; // multiply with ddij / di
					dSASA_4_neigh_dc2.y += sjk * dAdd * dddy;
					dSASA_4_neigh_dc2.z += sjk * dAdd * dddz;

					// add up component of the third term of the SASA itself
					Si3_t += sjk;
				}

				dSASA_4_neigh_dc.x = sji * dSASA_3_neigh_dc.x + dSASA_4_neigh_dc2.x;
				dSASA_4_neigh_dc.y = sji * dSASA_3_neigh_dc.y + dSASA_4_neigh_dc2.y;
				dSASA_4_neigh_dc.z = sji * dSASA_3_neigh_dc.z + dSASA_4_neigh_dc2.z;

				dSASA_3_neigh_dc.x += dSASA_3_neigh_dc2.x;
				dSASA_3_neigh_dc.y += dSASA_3_neigh_dc2.y;
				dSASA_3_neigh_dc.z += dSASA_3_neigh_dc2.z;

				// use the obtained factor to get first half of fourth term of derivative
				// and add it to the second half, already calculated in the first k loop and
				// already stored in dSASA_4_neigh_dc.
				dSASA_4_neigh_dc.x += dSASA_2_neigh_dc.x * Si3ik_t;
				dSASA_4_neigh_dc.y += dSASA_2_neigh_dc.y * Si3ik_t;
				dSASA_4_neigh_dc.z += dSASA_2_neigh_dc.z * Si3ik_t;

				// multiply all the derivative terms by their respective parameters (P2-P4)
				// of atom j as well as the solvatio parameter of the atom j
				dSASA_2_neigh_dc.mul(SASAatom[ja].P2 * SASAatom[ja].sigma); // multiply third term in change in SASA/dxi by its factor
				dSASA_3_neigh_dc.mul(SASAatom[ja].P3 * SASAatom[ja].sigma); // multiply third term in change in SASA/dxi by its factor
				dSASA_4_neigh_dc.mul(SASAatom[ja].P4 * SASAatom[ja].sigma); // multiply third term in change in SASA/dxi by its factor

				// and add them to the grand derivative sum of atom i
				SASAatom[i].SASAderiv.add(dSASA_2_neigh_dc);
				SASAatom[i].SASAderiv.add(dSASA_3_neigh_dc);
				SASAatom[i].SASAderiv.add(dSASA_4_neigh_dc);

				// add the SASA components together
				Si4 += (Si3_t * Si2_t);
				Si3 += Si3_t;
				Si2 += Si2_t;

				// calculate the self terms of the derivative of atom i (i.e. dAi/dxi)
				dAijdc_2.add(dAijdc_2t);
				dAijdc_4.x += Si3_t * dAijdc_2t.x;
				dAijdc_4.y += Si3_t * dAijdc_2t.y;
				dAijdc_4.z += Si3_t * dAijdc_2t.z;
			}

			SAi = SASAatom[i].P1 * Si + SASAatom[i].P2 * Si2 + SASAatom[i].P3 * Si3 + SASAatom[i].P4 * Si4;


			// multiply the self terms of derivative by their parameters and solvation parameters
			dAijdc_2.mul(SASAatom[i].P2 * SASAatom[i].sigma);
			dAijdc_4.mul(SASAatom[i].P4 * SASAatom[i].sigma);

			// dSASAi/dxi = SdAijdc_2 + SdAijdc_4;
			SASAatom[i].SASAderiv.add(dAijdc_2);
			SASAatom[i].SASAderiv.add(dAijdc_4);

			if(SAi < 0) { // LCPO can yield negative SASA for highly buried atoms, so prevent negative surface areas and forces
				SAi = 0;
				SASAatom[i].SASAderiv.zero();
			}

			// save individul atom SASAs
			SASAatom[i].SASA = SAi;
			// add things to wspace
			epot_cav += (double) SASAatom[i].SASA * SASAatom[i].sigma;
			wspace.atom[i].epot += (double) SASAatom[i].SASA * SASAatom[i].sigma;
			forcem.setTo(SASAatom[i].SASAderiv);
			forcem.mul(1 / (PhysicsConst::Angstrom));
			//printf("Analytical: %3d %e %e %e \n", i, forcem.x,forcem.y,forcem.z );

			wspace.cur.atom[i].f.add(forcem);
			// sum up total SASA
			totalSASA += SAi;
		}

		wspace.ene.epot += epot_cav;
		return 0;
	}
} // namespace Physics



