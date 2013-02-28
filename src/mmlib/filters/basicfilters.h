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

#ifndef __BASIC_FILTERS
#define __BASIC_FILTERS

#include <vector>

#include "forcefields/ffbonded.h" // Defines the 'Physics::Bond' class member
#include "tools/cloneholder.h" // Defines the 'CloneHolder' class member
#include "filterbase.h" // Defines the 'FilterBase' and 'MolFilterBase' base classes
#include "workspace/pospointer.h" // Defines the 'SelectedIndices' class memeber
#include "workspace/workspace.fwd.h"
#include "system/proximitygrid.h"

class MoleculeBase;
class BondOrder;

namespace Physics
{
	class Bond;
}

namespace Library
{
	class AngleSet;
}

//-------------------------------------------------
/// \brief  BRIEF DESCRIPTION
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
/// \author Jon Rea 
/// \todo STATE OF DEVELOPMENT
/// \bug BUGS?
class PD_API BondFilter : public FilterBase
{
public:
	BondFilter();
	BondFilter( const WorkSpace& _mol );

	virtual ~BondFilter(){}
	virtual BondFilter* clone() const { return new BondFilter(*this); }

	void setWSpace( const WorkSpace& _mol );
	void setPicker( const PickBase& _picker );
	void setDefaultPicker();

	void setTollerance( double d ); ///< Set the number of Angstroms over than the equilibrium distance a bond is allowed to be
	void setDefaultTollerance(); ///< Reassign the default tollerance level

	virtual bool passes(); ///< Does the current wspace look ok in terms of bond lengths?
	virtual std::string reason(); ///< What was the reason for the previous failure?

private:
	virtual void setInternalName();

	bool m_Setup;
	void setup();

	double m_Tollerance; ///< The number of Angstroms over than the equilibrium distance a bond is allowed to be

	const WorkSpace* m_WSpace;
	std::vector<Physics::Bond> m_Bonds;
	CloneHolder<PickBase> m_Picker;
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
class PD_API SingleBondFilter : public MolFilterBase
{
public:
	SingleBondFilter() 
		: MolFilterBase(),
		index(-1,-1,0.0,0.0)
	{
		name = "SingleBondFilter";
	}
	virtual ~SingleBondFilter(){}
	virtual SingleBondFilter* clone() const { return new SingleBondFilter(*this); }

	void setBond( const MoleculeBase& _mol, int i, int j, double _dIdeal, double _dTollerance );
	void setBond( const Physics::FF_Bonded& _BondedFF, int i, int j );
	void setBond( const Physics::FF_Bonded& _BondedFF, int i, int j, double _dTollerance );

	void setTollerance( double _Tollerance );

	virtual bool passes();
	virtual std::string reason();

protected:
	virtual void setInternalName();
	double forceConstantToTollerance( double _ForceConstant ) const;
	void basicAsserts() const;
	Physics::Bond index;
};


class PD_API AtomicSeparationFilter : public MolFilterBase
{
public:
	AtomicSeparationFilter() 
		: MolFilterBase(),
		index(-1,-1,0.0,0.0)
	{
		name = "AtomicSeparationFilter";
	}
	virtual ~AtomicSeparationFilter(){}
	virtual AtomicSeparationFilter* clone() const { return new AtomicSeparationFilter(*this); }

	void setBond( const MoleculeBase& _mol, int i, int j, double _dIdeal, double _dTollerance );

	void setTollerance( double _Tollerance );

	virtual bool passes();
	virtual std::string reason();

protected:
	void basicAsserts() const;
	virtual void setInternalName();
	Physics::Bond index;
};


//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Jon Rea 
///
/// \todo This filter is incomplete, requiring major work....
///
/// \bug BUGS?
///
class PD_API OmegaGroupFilter : public FilterBase
{
public:
	OmegaGroupFilter() 
		: m_SingleBondFilter()
	{
		AssessPrimaryDistance = true;
		AssessSecondaryDistances = true;
		AssessTorsion = true;
	}
	virtual ~OmegaGroupFilter(){}
	virtual OmegaGroupFilter* clone() const { return new OmegaGroupFilter(*this); }

	static const double peptideDist_CN;
	static const double peptideDist_ON;
	static const double peptideDist_CCa;
	static const double peptideDist_CaN;
	static const double peptideDist_CaCa;
	static const double peptideDist_OH;
	static const double peptideDist_CaH;
	static const double peptideDist_CH;

	void setTo( WorkSpace& wspace, size_t resIndex, double _tollerance = -1.0 );
	virtual bool passes();
	virtual std::string reason();

	bool AssessPrimaryDistance;
	bool AssessSecondaryDistances;
	bool AssessTorsion;

private:

	struct OmegaTorsion : public IndexQuartet
	{
		Maths::dvector* pI;
		Maths::dvector* pJ;
		Maths::dvector* pA;
		Maths::dvector* pB;
	};

	virtual void setInternalName();

	void addDistRestraint( const WorkSpace& wspace, int i, int j, double dIdeal );

	OmegaTorsion m_Torsion;
	SingleBondFilter m_SingleBondFilter;
	std::vector<AtomicSeparationFilter> m_SecondaryDistFilters;
};

class PD_API AngleSetFilter : public FilterBase
{
public:
	AngleSetFilter(){};
	virtual ~AngleSetFilter(){};
	virtual AngleSetFilter* clone() const { return new AngleSetFilter(*this); }

	void setTo( const MoleculeBase& mol, const Library::AngleSet& angleSet, const PickResidueBase& picker = PickEverything() );
	virtual bool passes();
	virtual std::string reason();

protected:
	const MoleculeBase* m_Mol;
	const Library::AngleSet* m_AngleSet;
	CloneHolder< PickResidueBase > m_Picker;
};

class PD_API ClashFilterBase : public FilterBase
{
public:
	ClashFilterBase();
	ClashFilterBase( MoleculeBase& _mol );
	ClashFilterBase( WorkSpace& _mol );
	virtual ~ClashFilterBase();

	void setOLapFac( double fac );

	void setMolecule( MoleculeBase& _mol );
	void setMolecule( WorkSpace& _mol );

	void setDefaultPickers(); ///< Optional call to reinitialise the pickers to their default values
	void setForPicker( const PickBase& _Picker );
	void setAgainstPicker( const PickBase& _Picker );

	virtual bool passes() = 0;
	virtual std::string reason() = 0;

protected:
	MoleculeBase* m_Mol; ///< A pointer to the underlying molecule
	const BondOrder* m_BondOrder;

	virtual void initCore() = 0; ///< Performs internal setup
	virtual void setInternalName();

	bool isProxyAtom( int i, int j ) const; ///< returns true if the atom pair are 1-2 1-3 or 1-4 bonded

	CloneHolder<PickBase> m_ForPicker; ///< Which atoms should we be testing for clashes?
	CloneHolder<PickBase> m_AgainstPicker; ///< Which atoms are we querying these against?

	SelectedIndices m_ForAtoms; ///< Cache the indexes of the 'for atoms' that we are querying the 'against atoms' against
	SelectedIndices m_AgainstAtoms; ///< Cache the indexes of the 'against atoms' that we are querying the 'for atoms' against

	/// The lowest allowable squared ratio of the square distance between two particles
	/// divided by the square sum of their VDW radii. Lower than this value is deemed
	/// to be severe a particle clash and therefore a filter failure.
	double m_SqrMinOverlapFactor; 

	/// Workspaces atoms cant be changed and therefore the index lookup only needs to be done once. 
	/// Molecule must calculate this list on every call to passes().
	bool m_AtomContainerIsVolatile; 
};

//-------------------------------------------------
/// \brief  Clash filter by picker
/// \details
/// Note: Use this class if your AGAINST atoms move lots. use ClashGridFilter if they dont.
/// The clash filter looks for atoms whos VDW radii overlap other atoms by more than a given ratio.
/// By default all heavy atoms are queried against all other heavy atoms, but this can be changed by using
/// different pickers for the clash test. 
/// \author Jon Rea 
class PD_API ClashFilter : public ClashFilterBase
{
public:
	ClashFilter();
	ClashFilter( MoleculeBase& _mol );
	ClashFilter( WorkSpace& _mol );
	virtual ~ClashFilter();
	virtual ClashFilter* clone() const { return new ClashFilter(*this); }

	virtual bool passes();
	virtual std::string reason();

protected:
	virtual void initCore(); ///< Performs internal setup

	virtual void setInternalName()
	{
		name = "ClashFilter";
	}
};

//-------------------------------------------------
/// \brief  Clash filter by proxy-grid
/// \details
/// Note: Use this class if your AGAINST atoms move little. Use ClashFilter if they do.
/// \author Jon Rea 
class PD_API ClashGridFilter : public ClashFilterBase
{
public:
	ClashGridFilter();
	ClashGridFilter( MoleculeBase& _mol );
	ClashGridFilter( WorkSpace& _mol );
	virtual ~ClashGridFilter();
	virtual ClashGridFilter* clone() const { return new ClashGridFilter(*this); }
	
	virtual bool passes();
	virtual std::string reason();

	double GridCutoff;

protected:
	virtual void initCore(); ///< Performs internal setup

	CloneHolder<ProximityGrid> m_ProxyGrid;

	virtual void setInternalName()
	{
		name = "ClashGridFilter";
	}
};

#endif

