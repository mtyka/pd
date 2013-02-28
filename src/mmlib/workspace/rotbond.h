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

#ifndef _ROTATABLE_BOND_H
#define _ROTATABLE_BOND_H

#include "componentbase.h"

#define MAXROTATABLEBONDSEGS 4

#include "workspace/workspace.fwd.h"
#include "system/molecule.fwd.h"


//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API RotatableBond{

public:
	enum RotBondType
	{
		Single,
		Planar,
		Double,
		Phi,
		Psi,
		Omega
	};

	RotatableBond(){
		i = -1;
		j = -1;
		for(int i = 0; i < MAXROTATABLEBONDSEGS; i++) {
			segmentstart[i] = -1;
			segmentend[i] = -1;
		}
		Type = Single;
		status = 0;
		status2 = 0;
	}
	int i, j; ///< two atoms on either side of the rotatable bond
	int ip, jp;///< two (arbitaryly chosen) atoms that make it into a torsion
	// these two atoms are used to calcalte *differences* in torsion angle only
	double lastphi;
	double phi;
	RotBondType Type;
	int status;
	int status2;
	int segmentstart[MAXROTATABLEBONDSEGS]; // segments of atoms to be rotated
	int segmentend[MAXROTATABLEBONDSEGS];
};


//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API RotBond: public WorkSpaceComponentBase
{
public:
	RotBond();

	virtual size_t memuse(int level);

	virtual int find( int atomI, int atomJ ); ///< returns the index of the rotatable bond between atoms i and j, or -1 if there isn't one.
	virtual int rotate(int irbond, double angle);
	virtual int rotate(int irbond, double angle, int maxresidue);
	virtual bool isDelocalised( int rotBondIndex ); 
	virtual double calcTorsion(int irbond);
	virtual void silence( Verbosity::Type _OutputLevel );

	virtual size_t size() const { return rotbond.size(); }
#ifndef SWIG
	virtual RotatableBond& operator[](size_t _index) { return rotbond[_index]; }
	virtual const RotatableBond& operator[](size_t _index) const { return rotbond[_index]; }
#endif

protected:
	virtual void reinit( WorkSpace* _wspace );
	virtual void clear();

private:
	std::vector<RotatableBond> rotbond;
	Verbosity::Type m_OutputLevel;
};


//-------------------------------------------------
//
/// \brief This is a dummy placeholder that does nothing and throws exceptions when 
/// used .. what a grumpy function in life .. :-)
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API RotBond_Dummy: public RotBond
{
public:
	RotBond_Dummy(){};

	virtual size_t memuse(int level){ return 0;};

	virtual int find( int atomI, int atomJ ){ throwError(); return 0; }
	virtual int rotate(int irbond, double angle){ throwError(); return 0; }
	virtual int rotate(int irbond, double angle, int maxresidue){ throwError(); return 0; }
	virtual bool isDelocalised( int rotBondIndex ){ throwError(); return 0; } 
	virtual double calcTorsion(int irbond){ throwError(); return 0; }
	virtual void silence( Verbosity::Type _OutputLevel ){}

	virtual size_t size() const { throwError(); return 0; }
#ifndef SWIG
	virtual RotatableBond& operator[](size_t _index) { throwError(); return dummy; } 
	virtual const RotatableBond& operator[](size_t _index) const { throwError(); return dummy; }
#endif

protected:
	virtual void reinit( WorkSpace* _wspace ){};
	virtual void clear(){};
	RotatableBond dummy;
	void throwError() const
	{
		throw(ProcedureException("Rotatable bond list service requested but no Rotatable bond list was loaded. use wspace.setRotatableBondList(..)" ) ); 
	}
};


//-------------------------------------------------
//
/// \brief Class describing a current rotateable bonds defining atoms and Scope of influence.
///
/// \details This class is used in Torsional Minimisation, SegBuilder and derived classes.
/// Unlike 'RotatableBond', this class acts upon an independent range of atoms, rather,
/// than an entire system. It will generate a chain break, and is used in specific
/// situaltions e.g. in the loop builder where the core atfer a given  section, and the core
/// before a given section must remain stationary,but the region between them can more.
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API RotationDefinition
{
public:
	MoleculeBase *molBase;

	Maths::dvector  *atom1; ///< anchor atom
	Maths::dvector  *atom2; ///< atom of the torsion bond
	Maths::dvector  *atom3; ///< atom of the torsion bond
	Maths::dvector  *atom4; ///< anchor atom

	// Atoms included between these indexes are rotated.
	int startAtomIndex; ///< index of the 1st atom to be effected by the rotation.
	int endAtomIndex; ///< index of the final atom to be effected by the rotation

	///< -- Minor HACK -- > (no real performance hit - we are only making one extra 'if' per rotation.
	///< When in 'reverse-rotation mode', the 'start to end' range does not cover the 'O' atom for the 
	///< omega rotation like it does for all other rotations - we have to make provisions to store its 
	///< index so that the roation actually moves the 'O' of the peptide group - grrrr!
	int arseAtom; 

	RotationDefinition():
	molBase( NULL ),
		startAtomIndex(-1),
		endAtomIndex(-1),
		arseAtom(-1),
		atom1( NULL ),
		atom2( NULL ),
		atom3( NULL ),
		atom4( NULL )
	{
	}

	inline double getCurrentTorsionAngle() { return calcTorsionAngle( *atom1, *atom2, *atom3, *atom4 ); }
	inline bool isValid() const { return startAtomIndex != -1; }

	void performRotation( double desiredAngle ); // sets the torsion to a given angle
	void perturbRotation( double deltaAngle ); // perturb the current torsion by a given angle change (delta)
};

#endif

