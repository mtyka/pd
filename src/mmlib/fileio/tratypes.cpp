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

#include "tratypes.h"

#include "forcefields/ffparam.h"
#include "forcefields/forcefield.h"
#include "workspace/workspace.h"

using namespace Physics;

namespace IO 
{
	BTF_Header::BTF_Header()
	{
		version = BTF_VERSION;

		time_t now = time(NULL);
		time(&now);
		tm* l_time = localtime(&now);
		dateyear = l_time->tm_year + 1900;
		datemonth = l_time->tm_mon + 1;
		dateday = l_time->tm_mday;

		strcpy(&ID0[0], "null");
		strcpy(&ID1[0], "null");
		strcpy(&ID2[0], "null");
		strcpy(&ID3[0], "null");

		char StandardAtomProperty[32][16] = {
			"Epsilon","Mass","Charge","",
			"","","","","","","","","","","","","","",
			"","","","","","","","","","","","","","" };
		memcpy(&customAtomProperty[0][0], &StandardAtomProperty[0][0], 32 * 16); // names of custom Properties

		char StandardEnergyEntry[64][16] = {
			"Step","time","cRMS","dRMS","etot","epot","ekin","epot_bond","epot_angle","epot_torsion",
			"epot_vdw","epot_elec","epot_surf","epot_pol","epot_pol_self","epot_pol_cross","epot_pol_tcross",
			"epot_pol_tot","epot_hb_total","epot_hb_native","epot_4","epot_5","epot_6","epot_custom",
			"data1", "data2", "data3", "data4", "data5", "data6", "data7", "data8", "data9", "data10",
			"", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
			"", "", "", "", "", "", "", "", "", "", "", "", "", "", "" };
		memcpy(&customEnergyEntry[0][0], &StandardEnergyEntry[0][0], 64 * 16); // names of custom energies

		memset(&descriptor[0], 0, 1024); // custom 1K descriptor ASCII area
		memset(&text[0], 0, 15 * 1024); // custom 15K text ASCII area
	}

	BTF_Header::~BTF_Header()
	{
		// Nothing to do
	}

	void BTF_Header::load(std::ifstream &_traFile)
	{
		_traFile.read((char*)this, sizeof(BTF_Header));
		if( version != BTF_VERSION )
		{
			THROW(ParseException,"The major version of the Tra_file and code_version do not match!");
		}
	}

	void BTF_Header::save(std::ofstream &_traFile, WorkSpace *_wspace, BFT_StandardIncludes _includeData, int _traStartOffset, int _additionalBlockSize)
	{
		Type = (int)(1 | _includeData); // 1 = AtomPos - these will always be defined!
		atoms = (int)_wspace->atom.size();
		residues = (int)_wspace->res.size();

		// determine blocksize
		blocksize = 4;
		blocksize += atoms * sizeof(float) * 3; // atom positions
		if( 0 != ( PhiPsis & _includeData) ) blocksize += residues * sizeof(float) * 2;
		if( 0 != ( Energies & _includeData) ) blocksize += sizeof(BTF_Energy);
		if( 0 != ( ForceVectors & _includeData) ) blocksize += atoms * sizeof(float) * 3;
		blocksize += _additionalBlockSize;

		trajectorystart = sizeof(BTF_Header) + 8 + (sizeof(BTF_SystemDefinitionEntry) * atoms) + _traStartOffset;

		// All done!
		_traFile.write((char*)this, sizeof(BTF_Header));
	}

	void BTF_Header::info( bool verbose )
	{
	}






	BTF_SystemDefinitionEntry::BTF_SystemDefinitionEntry()
	{
	}

	void BTF_SystemDefinitionEntry::info(bool verbose)
	{
		if(verbose)
		{
			printf("%4d %5s %5s %5s %5d %6s %8.1lf %8.1lf %8.1lf %1d %6.3lf %6.3lf\n",
				atomnumber, // atom number
				&pdbname[0], // pdb name
				&primitivetype[0], // primitive forcefield Type
				&altname[0], // alternative name
				parentnumber, // parent (i.e. residue) number
				&parentname[0], // parentname
				targetx, targety, targetz, // Xray strucuture coordinate if known
				structureknown, // 1 if valid atomic experimental data, 1 if unknown
				charge,
				radius );
		}
		else
		{
			printf("%4d %5s %5d %6s \n",
				atomnumber, // atom number
				&pdbname[0], // pdb name
				parentnumber, // parent (i.e. residue) number
				&parentname[0]); // parentname
		}
	}

	void BTF_SystemDefinitionEntry::setFrom(int _index, Particle* _particle)
	{
		// Atom
		atomnumber = _index;
		strncpy(&altname[0],_particle->rawname.c_str(), 8);
		strncpy(&pdbname[0],_particle->pdbname.c_str(), 8);
		strncpy(&primitivetype[0],_particle->type_name.c_str(), 8);

		// Residue
		strcpy(&parentname[0],_particle->parentname.c_str());
		parentnumber =_particle->ir; // parent (i.e. residue) number

		// Target
		targetx = (float)_particle->posRef().x; // Xray strucuture coordinate if known
		targety = (float)_particle->posRef().y;
		targetz = (float)_particle->posRef().z;
		structureknown = _particle->isKnownStructure() ? 1 : 0; // 1 if valid atomic experimental data, 1 if unknown

		// Atom Properties
		charge = (float)_particle->charge;
		radius = (float)_particle->radius;
		customProperty[0] = (float)_particle->epsilon;
		customProperty[1] = (float)_particle->mass;
		customProperty[2] = (float)_particle->charge;
		for(int n = 3; n < 32; n++)
		{
			customProperty[n] = 0.0; // currently not used
		}

		n_cov12atoms = (int)_particle->cov12atom.size();

		if(n_cov12atoms > 6)
			THROW(CodeException,"Trajectory::writeSystemDefinition() - Trajectory can hold a maximum of 6 covalent links per atom");

		for(int n = 0; n <_particle->cov12atom.size(); n++)
		{
			cov12atom[n] =_particle->cov12atom[n].i;
		}
	}

	void BTF_SystemDefinitionEntry::save(std::ofstream &_traFile)
	{
		_traFile.write((char*)this,sizeof(BTF_SystemDefinitionEntry));
	}

	void BTF_SystemDefinitionEntry::load(std::ifstream &_traFile)
	{
		_traFile.read((char*)this,sizeof(BTF_SystemDefinitionEntry));
	}






	BTF_Energy::BTF_Energy()
	{
		// blank these filler members...
		for(int i = 0; i < 30; i++)
		{
			fill[i] = 0.0f;
		}
	}

	void BTF_Energy::info(bool verbose)
	{
		printf("BTF_Energy::info() is not implemented (waiting for reform of WorkSpace energy storage...)\n");
	}

	void BTF_Energy::save(std::ofstream &_traFile, WorkSpace* _wspace)
	{
		Step = (float) _wspace->Step;
		time = (float) _wspace->Step; //1.0E15*(double)_wspace->param.Timestep * (double)_wspace->Step; // in seconds !
		cRMS = (float) _wspace->ene.cRMS;
		dRMS = (float) _wspace->ene.dRMS;

		// energy of current conformation
		etot = (float)(PhysicsConst::J2kcal * PhysicsConst::Na * _wspace->ene.etot);
		epot = (float)(PhysicsConst::J2kcal * PhysicsConst::Na * _wspace->ene.epot);
		ekin = (float)(PhysicsConst::J2kcal * PhysicsConst::Na * _wspace->ene.ekin);
		epot_bond = (float)(PhysicsConst::J2kcal * PhysicsConst::Na * _wspace->ene.epot_bond);
		epot_angle = (float)(PhysicsConst::J2kcal * PhysicsConst::Na * _wspace->ene.epot_angle);
		epot_torsion = (float)(PhysicsConst::J2kcal * PhysicsConst::Na * _wspace->ene.epot_torsion);
		epot_vdw = (float)(PhysicsConst::J2kcal * PhysicsConst::Na * _wspace->ene.epot_vdw);
		epot_elec = (float)(PhysicsConst::J2kcal * PhysicsConst::Na * _wspace->ene.epot_elec);
		epot_surf = (float)(PhysicsConst::J2kcal * PhysicsConst::Na * _wspace->ene.epot_surf);
		epot_pol = (float)(PhysicsConst::J2kcal * PhysicsConst::Na * _wspace->ene.epot_pol);
		epot_pol_self = (float)(PhysicsConst::J2kcal * PhysicsConst::Na * _wspace->ene.epot_pol_self);
		epot_pol_cross = (float)(PhysicsConst::J2kcal * PhysicsConst::Na * _wspace->ene.epot_pol_cross);
		epot_pol_totcross = (float)(PhysicsConst::J2kcal * PhysicsConst::Na * (_wspace->ene.epot_pol_cross + _wspace->ene.epot_elec));
		epot_pol_tot = (float)(PhysicsConst::J2kcal * PhysicsConst::Na * (_wspace->ene.epot_pol + _wspace->ene.epot_elec));
		epot_hb_total = 0;
		epot_hb_native = 0;
		epot_4 = 0;
		epot_5 = 0;
		epot_6 = 0;
		epot_custom = 0;

		data1 = (float)_wspace->ene.data1;
		data2 = (float)_wspace->ene.data2;
		data3 = (float)_wspace->ene.data3;
		data4 = (float)_wspace->ene.data4;
		data5 = (float)_wspace->ene.data5;
		data6 = (float)_wspace->ene.data6;
		data7 = (float)_wspace->ene.data7;
		data8 = (float)_wspace->ene.data8;
		data9 = (float)_wspace->ene.data9;
		data10 = (float)_wspace->ene.data10;

		_traFile.write((char*)this,sizeof(BTF_Energy));
	}

	void BTF_Energy::load(std::ifstream &_traFile)
	{
	}
}

