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

#ifndef __BTF_TYPES_H
#define __BTF_TYPES_H

#include <fstream>

#include "workspace/workspace.fwd.h"
class PD_API Particle;

namespace IO 
{

	enum BFT_StandardIncludes // Flagged enumeration
	{
		//AtomPos = 1 - Now a required statement, not a flag.
		PhiPsis = 2,
		// Rotamers = 4, Was here but is now deprecated...
		Energies = 8,
		ForceVectors = 16,

		MinimalIncludes = Energies,
		DefaultIncludes = MinimalIncludes,
		AllIncludes = PhiPsis | Energies | ForceVectors
	};






//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API BTF_Header
	{
	public:
		BTF_Header();
		~BTF_Header();

		void info(bool verbose);
		void load(std::ifstream &_traFile);
		void save(std::ofstream &_traFile, WorkSpace *_wspace, BFT_StandardIncludes _includeData, int _traStartOffset, int _additionalBlockSize);


		/// the internal file-format version
		static const int BTF_VERSION = 2; 


		///  version - currently 2. Used to tell any program parsing the file which specification version the file adheres to.
		int version; 

		///  flagged Type descriptor
		int Type; 

		///  nr of residues
		int residues; 

		///  nr of atoms
		int atoms; 


		// Blocksize of a trajectory entry. The size in bytes of the trajectory entry structure defined below. Note that although in the specification for any program that parses the file is a constant size, other padding data may be included by a client program per time step, which may have to be included in this region. The blocksize is therefore required to move to a particular position in the file. The number of steps that are stored in the file is not defined in the header. This can however be calculated from the blocksize and the filesize.
		int blocksize;

		///  start of first trajectory entry.  Again, as a client program may want to add custom padding data after the header, this field is required for any parsing program to know where the trajectory entries begin in the file.
		int trajectorystart; 
		int dateday;
		int datemonth;
		int dateyear;
		char ID0[32]; // custom ID strings
		char ID1[32];
		char ID2[32];
		char ID3[32];
		
		/// custom atom property strings – In the next section, all atoms are defined, including position properties and bonding (see below for the format). These will contain a number of custom floats (singles in VB, 4 bytes) depending on how the trajectory entry has been produced. For example a custom couloumbic radius would be defined here. Some specific names may be fixed here in the definition at a later point to allow any viewing programs to display properties in a context specific way.
		char customAtomProperty[32][16]; // names of custom Properties
		
		
		/// if energies are defined in the tra file per time step, their names must be defined here. All defined floats must be named. If these names are blank, even if a float is defined below in an energy structure, it will be ignored. 
		char customEnergyEntry[64][16]; 
		
		
		/// this will contain pre defined properties in the format.
		/// Property name = Property description
		/// e.g.
		/// ForceField=Amber94
		/// This will allow a folder full of tra files to be parsed for a particular property, e.g. where any files created between these dates with the following search method and property X enabled?
		char descriptor[1024]; 
		
		
		/// This 15k region will be human readable text. This will be displayed in a viewer, but is not to be machine parsed. It is purely to archive what method was used for referral by a user at a later date.
		char text[15 * 1024]; 
	};






//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API BTF_SystemDefinitionEntry
	{
	public:
		BTF_SystemDefinitionEntry();

		void load(std::ifstream &_traFile);
		void save(std::ofstream &_traFile);
		void setFrom(int _index, Particle* _particle);
		void info(bool verbose = false);

		/// The number of the particular atom. This doesn’t have to start from 0, be in order or be unique, however these would all be sensible. The atom number is internally handled by the viewer based on its position in the array. The integers for bonding (see below) refer to the position in the array, NOT to the atom number defined here.
		int atomnumber; 

		/// The standard name for a given atom that is defined by the PDB file format.
		char pdbname[8]; 

		/// The name that a given atom is given by a type based forcefield defining its properties.
		char primitivetype[8]; 

		/// Alternative name that is used internally by client programs – this property is parsed but disregarded by the viewer.
		char altname[8]; 


		/// parent (i.e. residue) number
		int parentnumber; 

		/// parentname, String of length 8 chars defining the name of the parent e.g. LYS – Lysine.
		char parentname[8]; 


		/// Xray strucuture coordinate if known
		float targetx, targety, targetz; 

		/// 1 if valid atomic experimental data, 1 if unknown
		int structureknown; 


		/// Defines atom 1 to atom 2 bonds from this atom. 4 integers defining the position in the array of the atoms to which this atom bonds. 
		int cov12atom[6]; 

		/// Defined the number of slots in the above array that are in use. 
		int n_cov12atoms; 

		/// Its is assumed that these are the two main properties that will always be defined, no matter what the forcefield, for a given protein. These are therefore explicitly named in the file format. 
		float charge; 

		///Its is assumed that these are the two main properties that will always be defined, no matter what the forcefield, for a given protein. These are therefore explicitly named in the file format.
		float radius; 

		/// These are extra per-atom properties that are defined in the file. The names of these were defined in the header. These floats are only parsed if a name has been defined for them in the header. These may be given some standardised names that the viewing application has intrinsic knowledge of in order to properly display the property to the user.
		float customProperty[32]; 
	};






//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API BTF_Energy
	{
	public:
		BTF_Energy();

		void info(bool verbose = false);
		void load(std::ifstream &_traFile);
		void save(std::ofstream &_traFile, WorkSpace* _wspace );

		float Step; // 0
		float time;
		float cRMS;
		float dRMS;

		float etot; // 4
		float epot;
		float ekin;
		float epot_bond;
		float epot_angle;
		float epot_torsion;
		float epot_vdw;
		float epot_elec;
		float epot_surf;
		float epot_pol;

		float epot_pol_self; // 14
		float epot_pol_cross;
		float epot_pol_totcross;
		float epot_pol_tot;
		float epot_hb_total;
		float epot_hb_native;
		float epot_4;
		float epot_5;
		float epot_6;
		float epot_custom;

		float data1;
		float data2;
		float data3;
		float data4;
		float data5;
		float data6;
		float data7;
		float data8;
		float data9;
		float data10;

		float fill[30]; //to make the structure's size 64*sizeof(float), reserving space for later additions
	};




} // namespace IO


#endif

