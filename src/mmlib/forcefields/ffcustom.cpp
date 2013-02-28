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

#include "forcefields/ffcustom.h"

#include "forcefields/ffparam.h"
#include "workspace/workspace.h"

using namespace Maths;
using namespace IO;

namespace Physics
{
	// -------------------------------------------------------------------------------
	// Class: CustomForce
	// -------------------------------------------------------------------------------
	// Individual forces in a user defined way - can be used for restraints etc..

	void CustomForce::info(const WorkSpace &wspace) const
	{ 
		// prints a little block of parameter information

		if(j >= 0) 
		{
			printf(" %6s(%3d):%4s %6s(%3d):%4s ",
				wspace.atom[i].parentl3name.c_str(),
				wspace.atom[i].ir,
				wspace.atom[i].pdbname.c_str(),
				wspace.atom[j].parentl3name.c_str(),
				wspace.atom[j].ir,
				wspace.atom[j].pdbname.c_str());
		} 
		else 
		{
			printf(" %6s(%3d):%4s (%6.1lf,%6.1lf,%6.1lf) ",
				wspace.atom[i].parentl3name.c_str(), 
				wspace.atom[i].ir, 
				wspace.atom[i].pdbname.c_str(), 
				absv.x,
				absv.y, 
				absv.z);
		}

		switch (Type) 
		{
		case HARMONIC: // dangerous - very strong forces at large separations
			printf(" Harmonic Tdist: %6.2lf Epsilon: %8.3lf \n", tdist, epsilon * PhysicsConst::J2kcal * PhysicsConst::Na);
			break;
		case BELLSHAPED:
			printf(" BellShaped Tdist: %6.2lf Epsilon: %8.3lf Gamma %4.1lf\n", tdist, epsilon * PhysicsConst::J2kcal * PhysicsConst::Na, gamma);
			break;
		case VSHAPED:
			printf(" VShaped Tdist: %6.2lf Epsilon: %8.3lf Gamma %4.1lf Beta %4.1lf\n",
				tdist, epsilon * PhysicsConst::J2kcal * PhysicsConst::Na, gamma, beta);
			break;
		case LINEARRESTRAINT:
			printf(" LinearRest. Tdist: %6.2lf Gamma: %8.3lf Tdist2: %6.2lf Delta: %8.3lf \n",
				tdist, gamma * PhysicsConst::J2kcal * PhysicsConst::Na, tdist2, delta * PhysicsConst::J2kcal * PhysicsConst::Na);
			break;
		default:
			break;
		}
	}

	void CustomForce::calcEnergyVerbose(FF_Custom * cff, ForcefieldBase::AtomicVerbosity level)
	{
		if(!Active)
			return;
	}

	void CustomForce::calcEnergy(FF_Custom * cff)
	{
		if(!Active)
			return;
	}

	void CustomForce::calcForce(FF_Custom * cff)
	{
		// Proxies
		WorkSpace& wspace = cff->getWSpace();

		double potential;
		double force;
		double totalforce = 0.0;
		double Dist_ij;
		double devx;
		dvector fv;

		if(!Active)
			return; // ignore inactive potentials

		if(j >= 0) 
		{
			fv.x =    double (wspace.cur.atom[i].p.x - wspace.cur.atom[j].p.x);
			fv.y =    double (wspace.cur.atom[i].p.y - wspace.cur.atom[j].p.y);
			fv.z =    double (wspace.cur.atom[i].p.z - wspace.cur.atom[j].p.z);
			//Dist_ij = (double) wspace.cur.atom[i].p.dist(wspace.cur.atom[j].p);
			Dist_ij = sqrt( sqr(fv.x) + sqr(fv.y) + sqr(fv.z) );
		} 
		else 
		{
			//Dist_ij = (double) wspace.cur.atom[i].p.dist(absv);
			fv.x = double (wspace.cur.atom[i].p.x - absv.x);
			fv.y = double (wspace.cur.atom[i].p.y - absv.y);
			fv.z = double (wspace.cur.atom[i].p.z - absv.z);
			Dist_ij = sqrt( sqr(fv.x) + sqr(fv.y) + sqr(fv.z) );
			if(Dist_ij < 0.001)
				Dist_ij = 0.001;
		}

		switch (Type) 
		{

		case HARMONIC:
			potential = 0.5 * epsilon * sqr(Dist_ij - tdist);
			force = epsilon * (Dist_ij - tdist) / PhysicsConst::Angstrom; //calculate force
			break;

		case BELLSHAPED:
			potential = epsilon / (1.0 + gamma * sqr(Dist_ij - tdist));
			force = -2.0 * epsilon * gamma * (Dist_ij - tdist) / sqr(1.0 + gamma * sqr(Dist_ij - tdist)) / PhysicsConst::Angstrom;
			// potential = epsilon / (1.0 + gamma * sqr(Dist_ij - tdist));
			// force = 2.0 * epsilon * gamma * (Dist_ij - tdist) / PhysicsConst::Angstrom / sqr(1.0 + gamma * sqr(Dist_ij - tdist));
			break;

		case VSHAPED:

			devx = beta * (Dist_ij - tdist) / gamma;

			if(devx > 8.0) 
			{
				potential = epsilon * gamma * devx;
				force = beta * epsilon / PhysicsConst::Angstrom;
			} 
			else if(devx < -8.0) 
			{
				potential = 0;
				force = 0;
			} 
			else 
			{
				potential = epsilon * gamma * log(1.0 + exp(devx));
				force = beta * (epsilon / (exp(-devx) + 1.0)) / PhysicsConst::Angstrom;
			}
			break;
		case LINEARRESTRAINT:
			potential = 0.0;
			force = 0.0;
			// left restraint (close side)
			if(Dist_ij < tdist) {
				potential += -gamma * (Dist_ij - tdist);
				force += -gamma / PhysicsConst::Angstrom; // constant repulsive force
			}
			// right restraint (far side)
			if(Dist_ij > tdist2) {
				potential += delta * (Dist_ij - tdist2);
				force += delta / PhysicsConst::Angstrom; // constant attarctive force
			}
		}

		wspace.ene.epot += potential;
		cff->epot += potential;
		// now apply the force

		fv.mul(force / Dist_ij);
		wspace.cur.atom[i].f.sub(fv);
		if(j >= 0)
			wspace.cur.atom[j].f.add(fv);
	}

	// this does the same calculation as above but returns the
	// pairwise energy
	void CustomForce::calcForceSingle(double Dist_ij, double &potential, double &force)
	{
		double devx;

		switch (Type) 
		{
		case HARMONIC:
			potential = 0.5 * epsilon * sqr(Dist_ij - tdist);
			force = epsilon * (Dist_ij - tdist) / PhysicsConst::Angstrom; //calculate force
			break;
		case BELLSHAPED:
			potential = epsilon / (1.0 + gamma * sqr(Dist_ij - tdist));
			force = -2.0 * epsilon * gamma * (Dist_ij - tdist) / sqr(1.0 + gamma * sqr(Dist_ij - tdist)) / PhysicsConst::Angstrom;
			break;

		case VSHAPED:
			devx = beta * (Dist_ij - tdist) / gamma;

			if(devx > 8.0) 
			{
				potential = epsilon * gamma * devx;
				force = beta * epsilon / PhysicsConst::Angstrom;
			} 
			else if(devx < -8.0) 
			{
				potential = 0;
				force = 0;
			} else 
			{
				potential = epsilon * gamma * log(1.0 + exp(devx));
				force = beta * (epsilon / (exp(-devx) + 1.0)) / PhysicsConst::Angstrom;
			}
			break;

		case LINEARRESTRAINT:
			potential = 0.0;
			force = 0.0;
			// left restraint (close side)
			if(Dist_ij < tdist) 
			{
				potential += -gamma * (Dist_ij - tdist);
				force += -gamma / PhysicsConst::Angstrom; // constant repulsive force
			}
			// right restraint (far side)
			if(Dist_ij > tdist2) 
			{
				potential += delta * (Dist_ij - tdist2);
				force += delta / PhysicsConst::Angstrom; // constant attarctive force
			}
		}
	}

	void CustomForce::printTextLine(WorkSpace &wspace, FILE * file)
	{
		if(j >= 0) 
		{
			fprintf(file, "%3d %4s %3d %4s ", wspace.atom[i].ir, wspace.atom[i].pdbname.c_str(), wspace.atom[j].ir, wspace.atom[j].pdbname.c_str());
		} 
		else 
		{
			fprintf(file, "%3d %4s (%.3lf,%.3lf,%.3lf) 0 ", wspace.atom[i].ir, wspace.atom[i].pdbname.c_str(), absv.x, absv.y, absv.z);
		}
		switch (Type) 
		{
		case HARMONIC:
			fprintf(file, "H %6.2lf %8.3lf \n", tdist, epsilon * PhysicsConst::J2kcal * PhysicsConst::Na);
			break;
		case BELLSHAPED:
			fprintf(file, "B %6.2lf %8.3lf %4.1lf\n", tdist, epsilon * PhysicsConst::J2kcal * PhysicsConst::Na, gamma);
			break;
		case VSHAPED:
			fprintf(file, "V %6.2lf %8.3lf %4.1lf %4.1lf\n", tdist, epsilon * PhysicsConst::J2kcal * PhysicsConst::Na, gamma, beta);
			break;
		case LINEARRESTRAINT:
			fprintf(file, "LR %6.2lf %8.3lf %6.2lf %8.3lf \n",
				tdist, gamma * PhysicsConst::J2kcal * PhysicsConst::Na, tdist2, delta * PhysicsConst::J2kcal * PhysicsConst::Na);
			break;
		}
	}

	bool CustomForce::compareAtoms(const CustomForce& cf)
	{
		if(j < 0) return false;
		if(cf.j < 0) return false;
		if((cf.i == i) && (cf.j == j))
			return true;
		if((cf.i == j) && (cf.j == i))
			return true;
		return false;
	}

	void FF_Custom::info() const
	{ 
		// Proxies
		const WorkSpace& wspace = getWSpace();

		// prints a little block of parameter information
		printf(" --FF_Custom -----------------\n");
		if(!Active)
			printf(" -----------------------<INACTIVE>---\n");

		printf(" TotalForces: %d \n\n", force.size());

		if(verbose != 0) 
		{
			for(int i = 0; i < force.size(); i++) 
			{
				force[i].info(wspace);
			}
		}
		printf("\n");
	}

	void FF_Custom::infoLine() const 
	{ 
		// prints a line of current energies
		printf("\t%6.1lf", double (epot) * PhysicsConst::J2kcal * PhysicsConst::Na);
	}

	void FF_Custom::infoLineHeader() const 
	{
		// prints the headers for the above function
		printf("\t%7s", name_infoheader.c_str());
	}


	void FF_Custom::calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level)
	{
		if(!Active)
			return;
		epot = 0.0;
		int nforces = force.size();
		calcForces();
		if(level == Summary) 
		{
			printf(" %-16s%10.3lf kcal/mol\n",
				name_verbose.c_str(),
				double (epot) * PhysicsConst::J2kcal * PhysicsConst::Na);
			return;
		}

		for(int i = 0; i < nforces; i++) 
		{
			force[i].calcEnergyVerbose(this, level);
		}
	}

	void FF_Custom::calcEnergies()
	{
		if(!Active)
			return;
		int nforces = force.size();
		epot = 0.0;
		for(int i = 0; i < nforces; i++) 
		{
			force[i].calcEnergy(this);
		}
	}

	void FF_Custom::calcForces()
	{
		if(!Active)
			return;
		int nforces = force.size();
		epot = 0.0;
		for(int i = 0; i < nforces; i++) 
		{
			force[i].calcForce(this);
		}
	}

	// calculates the force & energy at a certain i<->j distnace
	// assumes all the component forces act on the same two particles
	// otherwise the answer is meaningless
	void FF_Custom::calcForceSingle(double Dist_ij, double &totpot, double &totforce)
	{
		double ipot,iforce;
		int nforces = force.size();
		epot = 0.0;
		totpot=0;
		totforce=0;
		for(int i = 0; i < nforces; i++) 
		{
			force[i].calcForceSingle(Dist_ij, ipot, iforce);
			totpot += ipot;
			totforce += iforce;
		}
	}

	// these (de)activate all or subsets of the custom forces (according to ID)
	void FF_Custom::activateforces()
	{
		for(int a = 0; a < force.size(); a++)
			force[a].activate();
	}
	void FF_Custom::activateforces(int ID)
	{
		for(int a = 0; a < force.size(); a++)
			if(force[a].ID == ID)
				force[a].activate();
	}
	void FF_Custom::activateforcesnot(int ID)
	{
		for(int a = 0; a < force.size(); a++)
			if(force[a].ID != ID)
				force[a].activate();
	}

	void FF_Custom::deactivateforces()
	{
		for(int a = 0; a < force.size(); a++)
			force[a].deactivate();
	}
	void FF_Custom::deactivateforces(int ID)
	{
		for(int a = 0; a < force.size(); a++)
			if(force[a].ID == ID)
				force[a].deactivate();
	}
	void FF_Custom::deactivateforcesnot(int ID)
	{
		for(int a = 0; a < force.size(); a++)
			if(force[a].ID != ID)
				force[a].deactivate();
	}

	int FF_Custom::findForce(const CustomForce& cf)
	{
		for(int a = 0; a < force.size(); a++)
			if(force[a].compareAtoms(cf))
				return a;
		return -1;
	}



	// loads&parses a file containing the restraint definitions
	int FF_Custom::loadTextFile(const char *filename)
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		char buffer[256];
		char command[256];
		char atomname1[30];
		char atomname2[30];
		int ati, atj;
		int residue1, residue2;
		char potentialtype[30];
		int beginning[6];
		int end[6];
		int nelements;
		int npot = 0;
		double param1, param2, param3, param4;

		bool insideListBlock = false;

		TextFilePtr file(filename, "r");
		if(!file.loaded()) 
		{
			printf("ERROR: File not found: %s \n", filename);
			return -1;
		}

		int line = 1;
		int lines;
		printf("reading custom force file .. %s \n", filename);
		while(!feof((FILE *) file)) 
		{
			if((lines = file.readline(&buffer[0], 256)) < 0)
				break;
			line += lines;

			sscanf(&buffer[0], "%s", &command[0]);

			if(strcmp(&command[0], "CUSTOMFORCELIST") == 0) {
				insideListBlock = true;
			}

			if(!insideListBlock)
				continue;
			if(strcmp(&command[0], "ENDCUSTOMFORCELIST") == 0) {
				insideListBlock = false;
				continue;
			}

			nelements = chopString(&buffer[0], &beginning[0], &end[0], 4);
			if(nelements < 4)
				continue;
			if(sscanf(&buffer[beginning[0]], "%d", &residue1) != 1)
				continue;
			if(sscanf(&buffer[beginning[2]], "%d", &residue2) != 1)
				continue;

			strncpy(&atomname1[0], &buffer[beginning[1] + 1], -beginning[1] + end[1]);
			atomname1[-beginning[1] + end[1]] = 0;
			strncpy(&atomname2[0], &buffer[beginning[3] + 1], -beginning[3] + end[3]);
			atomname2[-beginning[3] + end[3]] = 0;

			if(verbose != 0) 
			{
				printf("read line: -%d- -%s- -%d- -%s-\n", residue1, &atomname1[0], residue2, &atomname2[0]);
			}
			// find the two particles of concern
			ati = wspace.findParticle(residue1, &atomname1[0]);
			if(ati < 0) 
			{
				printf("ERROR: Particle 1 not found, %s, line %d\n", filename, line);
				return -1;
			}

			CustomForce *cf;

			if(atomname2[0] != '(') 
			{
				atj = wspace.findParticle(residue2, &atomname2[0]);
				if(atj < 0) {
					printf("ERROR: Particle 2 not found, %s, line %d\n", filename, line);
					return -1;
				}
				if(ati == atj) { // force cant act between two identical atoms
					printf("ERROR: Particle 1 & 2 are identical, %s, line %d\n", filename, line);
					return -1;
				}
				cf = new CustomForce(ati, atj);
			} 
			else 
			{
				dvector absv;

				for(int c = 0; c < (int) strlen(&atomname2[0]); c++) 
				{
					if(atomname2[c] == '(')
						atomname2[c] = ' ';
					if(atomname2[c] == ',')
						atomname2[c] = ' ';
					if(atomname2[c] == ')')
						atomname2[c] = ' ';
				}

				double x, y, z;
				sscanf(&atomname2[0], "%lf %lf %lf", &x, &y, &z);
				absv.setTo(x, y, z);

				cf = new CustomForce(ati, absv);
			}
			// Now create custom force

			if(sscanf(&buffer[end[3] + 1], "%s", &potentialtype[0]) != 1)
			{
				printf("Potential Type Expected, %s, line %d\n", filename, line);
				return -1;
			}

			if(strcmp(&potentialtype[0], "H") == 0) 
			{
				if(sscanf(&buffer[end[3] + 1], "%s %lf %lf", &potentialtype[0], &param1, &param2) < 3) 
				{
					printf("2 parameters expected after %s, %s, line %d\n", &potentialtype[0], filename, line);
					return -1;
				}
				cf->setToHarmonic(param1, param2 * PhysicsConst::kcal2J / PhysicsConst::Na);
			} 
			else if(strcmp(&potentialtype[0], "B") == 0) 
			{
				if(sscanf(&buffer[end[3] + 1], "%s %lf %lf %lf", &potentialtype[0], &param1, &param2, &param3) < 4) {
					printf("3 parameters expected after %s, %s, line %d\n", &potentialtype[0], filename, line);
					return -1;
				}
				cf->setToBellShaped(param1, param2 * PhysicsConst::kcal2J / PhysicsConst::Na, param3);
			} 
			else if(strcmp(&potentialtype[0], "V") == 0) {
				if(sscanf(&buffer[end[3] + 1], "%s %lf %lf %lf %lf", &potentialtype[0], &param1, &param2, &param3, &param4) < 5) {
					printf("4 parameters expected after %s, %s, line %d\n", &potentialtype[0], filename, line);
					return -1;
				}
				cf->setToVShaped(param1, param2 * PhysicsConst::kcal2J / PhysicsConst::Na, param3, param4);
			} 
			else if(strcmp(&potentialtype[0], "LR") == 0) 
			{
				if(sscanf(&buffer[end[3] + 1], "%s %lf %lf %lf %lf", &potentialtype[0], &param1, &param2, &param3, &param4) < 5) 
				{
					printf("4 parameters expected after %s, %s, line %d\n", &potentialtype[0], filename, line);
					return -1;
				}
				cf->setToLinearRestraint(param1, param2 * PhysicsConst::kcal2J / PhysicsConst::Na, param3, param4 * PhysicsConst::kcal2J / PhysicsConst::Na);
			} 
			else 
			{
				printf("Potential Type %s unknown, %s, line %d\n", &potentialtype[0], filename, line);
				return -1;
			};

			addForce(*cf);
			delete cf;
			npot++;
		}

		printf("read %d potentials \n", npot);
		return 0;
	}

	// loads&parses a file containing the restraint definitions
	int FF_Custom::saveTextFile(const char *filename)
	{
		// Proxies
		WorkSpace& wspace = getWSpace();

		TextFilePtr file(filename, "w");
		fprintf(file, "CUSTOMFORCELIST\n");
		for(int a = 0; a < force.size(); a++)
			force[a].printTextLine(wspace,file);
		fprintf(file, "ENDCUSTOMFORCELIST\n");
		return -1;
	}

	int FF_Custom::load(const char *filename)
	{ 
		// loads the forcefield from a binary form
		FILE *file;
		int n;
		int i;
		CustomForce cf;

		file = fopen(filename, "rb");
		if(file == NULL)
			return -1;

		// write number of individual forces
		fread((void *) &n, sizeof(int), 1, file);

		for(i = 0; i < n; i++) {
			fread((void *) &cf, sizeof(CustomForce), 1, file);
			force.push_back(cf);
		}

		if(force.size() != n) {
			printf("ERROR: Error while loading CustomForce file\n");
		}

		fclose(file);
		return 0;
	}


	int FF_Custom::save(const char *filename)
	{ 
		// saves the forcefield in a binary form
		FILE *file;
		int n = force.size();
		int i;

		file = fopen(filename, "wb");

		// write number of individual forces
		fwrite((void *) &n, sizeof(int), 1, file);

		for(i = 0; i < n; i++) 
		{
			fwrite((void *) &force[i], sizeof(CustomForce), 1, file);
		}

		fclose(file);
		return 0;
	}
}

