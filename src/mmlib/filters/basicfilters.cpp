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
#include "maths/maths.h"
#include "system/molecule.h"
#include "workspace/workspace.h"
#include "workspace/bondorder.h"
#include "pickers/basicpickers.h"
#include "library/angleset.h"
#include "basicfilters.h"

BondFilter::BondFilter() : 
	FilterBase()
{
	setDefaultTollerance();
	setDefaultPicker();
	m_Setup = false;
	m_WSpace = NULL;
}

BondFilter::BondFilter( const WorkSpace& _wspace ) :
	FilterBase()
{
	setDefaultTollerance();
	setDefaultPicker();
	setWSpace( _wspace );
}

void BondFilter::setInternalName()
{
	name = "BondFilter";
}

void BondFilter::setWSpace( const WorkSpace& _wspace )
{
	m_WSpace = &_wspace;
	m_Setup = false;
}

void BondFilter::setDefaultPicker()
{
	m_Setup = false;
	m_Picker = PickAllAtoms();
}

void BondFilter::setPicker( const PickBase& _picker )
{
	m_Setup = false;
	m_Picker = _picker;
}

void BondFilter::setTollerance( double d ) 
{ 
	m_Setup = false;
	m_Tollerance = d; 
}

void BondFilter::setDefaultTollerance() 
{ 
	m_Setup = false;
	m_Tollerance = 0.5; // Angstroms more than the equilibrium distance
}

void BondFilter::setup()
{
	if( m_WSpace == NULL ) throw NullArgumentException("BondFilter internal wspace pointer is NULL. Did you forget to call setWSpace()?");
	if( !m_Picker.assigned() ) throw NullArgumentException("BondFilter internal picker pointer is NULL");

	// Set up from the bond-information contained within the wspace
	// Proxies
	const WorkSpace& wspace = *m_WSpace;
	const ParticleStore& atomparam = wspace.atom;

	m_Bonds.clear(); //set counters to zero

	for(size_t i = 0; i < wspace.atom.size(); i++) 
	{
		if( !m_Picker->matches( wspace.atom[i] ) ) continue;

		//collect all the covalent links, assume they're all bond...
		for(int c = 0; c < atomparam[i].cov12atom.size(); c++) 
		{ 
			int j = atomparam[i].cov12atom[c].i;
			if(i < j) continue; // only record half (i.e. don't doubly record 1-3 and 3-1)
			if( !m_Picker->matches( wspace.atom[j] ) ) continue;

			Physics::Bond newbond;

			newbond.i = (int)i; //record atoms
			newbond.j = j;

			int bindex = wspace.ffps().findBondType(atomparam[newbond.i].FFType, atomparam[newbond.j].FFType);

			ASSERT( bindex >= 0, ProcedureException, "BondFilter:: Could not find ff bond type");

			// An angletype was found successfully, load parameters
			newbond.k = wspace.ffps().BondType[bindex].forceconstant; // force constant
			newbond.l = wspace.ffps().BondType[bindex].length; // eq. bond length
			newbond.l += m_Tollerance; // add our tollerance level to the equlibrium distance
			newbond.l *= newbond.l; // Make it internally into a *SQUARE* distance

			m_Bonds.push_back(newbond);
		}
	}

	m_Setup = true;
}

bool BondFilter::passes()
{
	if( !m_Setup ) setup();

	const WorkSpace& wspace = *m_WSpace;
	const SnapShotAtom *atom = wspace.cur.atom;

	for( size_t i = 0; i < m_Bonds.size(); i++ )
	{
		Physics::Bond& bond = m_Bonds[i];
		double sqrdist = atom[bond.i].p.sqrdist(atom[bond.j].p);
		if( sqrdist > bond.l ) 
			return false;
	}

	return true;
}

std::string BondFilter::reason()
{
	if( !m_Setup ) setup();

	const WorkSpace& wspace = *m_WSpace;
	const SnapShotAtom *atom = wspace.cur.atom;

	for( size_t i = 0; i < m_Bonds.size(); i++ )
	{
		Physics::Bond& bond = m_Bonds[i];
		double sqrdist = atom[bond.i].p.sqrdist(atom[bond.j].p);
		if( sqrdist > bond.l ) 
		{
			StringBuilder sb( 80 );
			sb.setFormat("Bond too long: Atom i: %d, Atom j %d, Tollerance Length: %8.3f, Actual Length: %8.3f")
				(bond.i)(bond.j)(sqrt(bond.l))(sqrt(sqrdist));
			return sb.toString();
		}
	}

	return "Pass";
}

void SingleBondFilter::setInternalName()
{
	name = "SingleBondFilter";
}

void SingleBondFilter::basicAsserts() const
{
	ASSERT( m_Mol != NULL, CodeException, "Internal NULL reference found in SingleBondFilter");
	ASSERT( index.i >= 0 && index.i < m_Mol->atom.size(), ArgumentException, "Bond index 'i' must be within the range of the molecule." );
	ASSERT( index.j >= 0 && index.j < m_Mol->atom.size(), ArgumentException, "Bond index 'j' must be within the range of the molecule." );
	ASSERT( index.l >= 0.0, ArgumentException, "SingleBondFilter ideal distance must be >0.0" );
	ASSERT( index.k >= 0.0, ArgumentException, "SingleBondFilter tollerance distance delta must be >0.0" );
}

void SingleBondFilter::setBond( const MoleculeBase& _mol, int i, int j, double _dIdeal, double _dTollerance )
{
	m_Mol = &_mol;
	index.i = i;
	index.j = j;
	index.l = _dIdeal;
	index.k = _dTollerance;
	basicAsserts(); // assert that we have just set something sensible
}

double SingleBondFilter::forceConstantToTollerance( double _ForceConstant ) const
{
	double tollerance = _ForceConstant / 3.5e-18; // very empiricle!
	return tollerance;
}

void SingleBondFilter::setBond( const Physics::FF_Bonded& _BondedFF, int i, int j )
{
	m_Mol = &_BondedFF.getWSpace();
	index = _BondedFF.getBond( i, j ); // this will throw an exception on failure
	index.k = forceConstantToTollerance( index.k );
	basicAsserts(); // assert that we have just set something sensible
}

void SingleBondFilter::setBond( const Physics::FF_Bonded& _BondedFF, int i, int j, double _dTollerance )
{
	setBond( _BondedFF, i, j );
	index.k = _dTollerance;
	ASSERT( index.k >= 0.0, ArgumentException, "SingleBondFilter tollerance distance delta must be >0.0" );
}

void SingleBondFilter::setTollerance( double _Tollerance )
{
	index.k = _Tollerance;
	ASSERT( index.k >= 0.0, ArgumentException, "SingleBondFilter tollerance distance delta must be >0.0" );
}

bool SingleBondFilter::passes()
{
	ASSERT( m_Mol != NULL, CodeException, "Internal NULL reference found in SingleBondFilter. Has setBond() been called?");
	double dist = m_Mol->calcAtomDistance( index.i, index.j );
	double absDiff = fabs( dist - index.l );
	return absDiff < index.k;
}

std::string SingleBondFilter::reason()
{
	ASSERT( m_Mol != NULL, CodeException, "Internal NULL reference found in SingleBondFilter. Has setBond() been called?");
	
	double dist = m_Mol->calcAtomDistance( index.i, index.j );
	double deviation = fabs( dist - index.l );

	StringBuilder sb;
	if( deviation < index.k )
	{
		sb.setFormat("Bond %d-%d is within tollerance. Length:%5.2lf Ideal:%5.2lf Deviation:%6.2lf Tollerance:%5.2lf")
			(index.i)(index.j)(dist)(index.l)(deviation)(index.k);
	}
	else
	{
		sb.setFormat("Bond %d-%d is NOT within tollerance. Length:%5.2lf Ideal:%5.2lf Deviation:%6.2lf Tollerance:%5.2lf")
			(index.i)(index.j)(dist)(index.l)(deviation)(index.k);
	}

	return sb.toString();
}

void AtomicSeparationFilter::setInternalName()
{
	name = "AtomicSeparationFilter";
}

void AtomicSeparationFilter::basicAsserts() const
{
	ASSERT( m_Mol != NULL, CodeException, "Internal NULL reference found in SingleBondFilter");
	ASSERT( index.i >= 0 && index.i < m_Mol->atom.size(), ArgumentException, "Bond index 'i' must be within the range of the molecule." );
	ASSERT( index.j >= 0 && index.j < m_Mol->atom.size(), ArgumentException, "Bond index 'j' must be within the range of the molecule." );
	ASSERT( index.l >= 0.0, ArgumentException, "SingleBondFilter ideal distance must be >0.0" );
	ASSERT( index.k >= 0.0, ArgumentException, "SingleBondFilter tollerance distance delta must be >0.0" );
}

void AtomicSeparationFilter::setBond( const MoleculeBase& _mol, int i, int j, double _dIdeal, double _dTollerance )
{
	m_Mol = &_mol;
	index.i = i;
	index.j = j;
	index.l = _dIdeal;
	index.k = _dTollerance;
	basicAsserts(); // assert that we have just set something sensible
}

void AtomicSeparationFilter::setTollerance( double _dTollerance )
{
	index.k = _dTollerance;
	basicAsserts(); // assert that we have just set something sensible
}

bool AtomicSeparationFilter::passes()
{
	return fabs( m_Mol->atomxyz( index.i ).dist( m_Mol->atomxyz( index.j ) ) - index.l ) < index.k;
}

std::string AtomicSeparationFilter::reason()
{
	double dist =  fabs( m_Mol->atomxyz( index.i ).dist( m_Mol->atomxyz( index.j ) ) );
	double deviation = dist - index.l;
	if( deviation < index.k )
	{
		return "Bond within tollerance";
	}
	else
	{
		StringBuilder sb;
		sb.setFormat("Separation %d-%d is NOT within tollerance. Length:%5.2lf Ideal:%5.2lf Deviation:%6.2lf Tollerance:%5.2lf")
			(index.i)(index.j)(dist)(index.l)(deviation)(index.k);
		return sb.toString();
	}
}

const double OmegaGroupFilter::peptideDist_CN = 1.33;
const double OmegaGroupFilter::peptideDist_ON = 2.25;
const double OmegaGroupFilter::peptideDist_CCa = 2.44;
const double OmegaGroupFilter::peptideDist_CaN = 2.43;
const double OmegaGroupFilter::peptideDist_CaCa = 3.80;
const double OmegaGroupFilter::peptideDist_OH = 3.02;
const double OmegaGroupFilter::peptideDist_CaH = 2.49;
const double OmegaGroupFilter::peptideDist_CH = 1.91;

void OmegaGroupFilter::setInternalName()
{
	name = "OmegaGroupFilter";
}

const double DEFAULT_OMEGA_DIST_TOLLERANCE = 0.4; // 0.4 seems an ok sort of test ... :-S
const double DEFAULT_OMEGA_TORSION_TOLLERANCE = 0.436; // 0.436radians == 25degrees ---> My data says this is absolute !MAX! for real PDB data

void OmegaGroupFilter::addDistRestraint( const WorkSpace& wspace, int i, int j, double dIdeal )
{
	if( i == -1 || j == -1 )
		return;
	AtomicSeparationFilter dist;
	dist.setBond( wspace, i, j, dIdeal, DEFAULT_OMEGA_DIST_TOLLERANCE ); 
	m_SecondaryDistFilters.push_back(dist);
}

void OmegaGroupFilter::setTo( WorkSpace& wspace, size_t resIndex, double _tollerance )
{
	if( wspace.isMolStartResIndex(resIndex) )
	{
		StringBuilder sb;
		sb.setFormat("You cannot ask for the omega angle of residue '%d'. As it is the start of a chain it doesnt have one by definition.")(resIndex);
		throw ArgumentException( sb.toString() );
	}
	Residue& res1 = wspace.res[resIndex-1];
	Residue& res2 = wspace.res[resIndex];
	int i = res1.iC;
	int j = res2.iN;
	ASSERT( i != -1, ProcedureException, "Could not find a valid omega 'C' atom");
	ASSERT( j != -1, ProcedureException, "Could not find a valid omega 'N' atom");
	m_SingleBondFilter.setBond( wspace, i, j, peptideDist_CN, _tollerance > 0.0 ? _tollerance : DEFAULT_OMEGA_DIST_TOLLERANCE );

	m_SecondaryDistFilters.clear();
	// addDistRestraint( wspace, res1.iC, res2.iN, peptideDist_CN ); - this is the m_SingleBondFilter
	addDistRestraint( wspace, res1.iO, res2.iN, peptideDist_ON );
	addDistRestraint( wspace, res1.iC, res2.iCA, peptideDist_CCa );
	addDistRestraint( wspace, res1.iCA, res2.iN, peptideDist_CaN );
	addDistRestraint( wspace, res1.iCA, res2.iCA, peptideDist_CaCa );
	addDistRestraint( wspace, res1.iO, res2.iH, peptideDist_OH );
	addDistRestraint( wspace, res1.iCA, res2.iH, peptideDist_CaH );
	addDistRestraint( wspace, res1.iC, res2.iH, peptideDist_CH );

	m_Torsion.a = res1.iCA;
	m_Torsion.i = res1.iC;
	m_Torsion.j = res2.iN;
	m_Torsion.b = res2.iCA;
	m_Torsion.pA = &wspace.atomxyz(m_Torsion.a);
	m_Torsion.pI = &wspace.atomxyz(m_Torsion.i);
	m_Torsion.pJ = &wspace.atomxyz(m_Torsion.j);
	m_Torsion.pB = &wspace.atomxyz(m_Torsion.b);
}

bool OmegaGroupFilter::passes()
{
	if( AssessPrimaryDistance )
	{
		if( !m_SingleBondFilter.passes() )
		{
			return false;
		}
	}
	if( AssessSecondaryDistances )
	{
		for( size_t i = 0; i < m_SecondaryDistFilters.size(); i++ )
		{
			if( !m_SecondaryDistFilters[i].passes() )
			{
				return false;
			}
		}
	}
	if( AssessTorsion )
	{
		const double halfwobble = DEFAULT_OMEGA_TORSION_TOLLERANCE; // plus minus 30.0 degrees is deemed ok
		double tor = Maths::calcTorsionAngle( *m_Torsion.pA, *m_Torsion.pI, *m_Torsion.pJ, *m_Torsion.pB );
		if( (Maths::MathConst::PI - halfwobble) <= tor && tor <= Maths::MathConst::PI )
		{
			// trans
			return true;
		}
		else if( (-Maths::MathConst::PI + halfwobble) >= tor && tor >= -Maths::MathConst::PI )
		{
			// trans
			return true;
		}
		else if( -halfwobble <= tor && tor <= halfwobble )
		{
			// cis
			return true;
		}
		else
		{
			return false; // badness
		}
	}
	return true;
}

std::string OmegaGroupFilter::reason()
{
	return "Unknown";
}

void AngleSetFilter::setTo( const MoleculeBase& mol, const Library::AngleSet& angleSet, const PickResidueBase& picker )
{
	m_Mol = &mol;
	m_AngleSet = &angleSet;
	m_Picker = picker;
}

bool AngleSetFilter::passes()
{
	const double DEFAULT_SEARCH_RADIUS = 50.0;
	const MoleculeBase& mol = *m_Mol;
	for( size_t i = 0; i < mol.nResidues(); i++ )
	{
		if( m_Picker->matches( mol.res[i] ) )
		{
			double phi;
			double psi;
			double omega;
			mol.calcResiduePhiPsi( i, phi, psi );
			mol.calcResidueOmega( i, omega );

			// 180-degree translation required!!!! Due to angleset torsion definitions!!
			psi = psi - Maths::MathConst::PI;
			Library::AngleSet::EnsureRadianRange( psi );
			omega = omega - Maths::MathConst::PI;
			Library::AngleSet::EnsureRadianRange( omega );

			double c_phi;
			double c_psi;
			double c_omega;
			m_AngleSet->getClosestAngleGroup( m_Mol->getSequence().getSingleResLetter(i),
				phi, c_phi, psi, c_psi, omega, c_omega );

			double omgDiff = c_omega - omega;
			if( omgDiff < 0.0 ) omgDiff = -omgDiff;
			if( omgDiff > Maths::MathConst::PI ) omgDiff = Maths::MathConst::TwoPI - omgDiff;
			if( omgDiff > DEFAULT_SEARCH_RADIUS ) // DEFAULT_OMEGA_TORSION_TOLLERANCE is to rigid for what i want this for, even though "correct", will think about it more later...
			{
				return false;
			}

			const double SqrSearchRadius = Maths::sqr( Maths::DegToRad( DEFAULT_SEARCH_RADIUS ) );
			
			// Ensure we take circular space into account!!
			double phiDiff = c_phi - phi;
			if( phiDiff < 0.0 ) phiDiff = -phiDiff;
			if( phiDiff > Maths::MathConst::PI ) phiDiff = Maths::MathConst::TwoPI - phiDiff;

			double psiDiff = c_psi - psi;
			if( psiDiff < 0.0 ) psiDiff = -psiDiff;
			if( psiDiff > Maths::MathConst::PI ) psiDiff = Maths::MathConst::TwoPI - psiDiff;

			if( SqrSearchRadius < ((phiDiff * phiDiff) + (psiDiff * psiDiff)) )
			{
				return false;
			}
		}
	}
	return true;
}

std::string AngleSetFilter::reason()
{
	return "Unknown";
}


const double DEFAULT_SQR_MAX_OVERLAP_FACTOR = std::pow(0.8,2);

ClashFilterBase::ClashFilterBase() 
	: FilterBase(), 
	m_Mol(NULL), 
	m_SqrMinOverlapFactor(DEFAULT_SQR_MAX_OVERLAP_FACTOR), 
	m_AtomContainerIsVolatile(true),
	m_BondOrder(NULL)
{
	setDefaultPickers();
}

ClashFilterBase::ClashFilterBase( MoleculeBase& _mol )
	: FilterBase(), 
	m_Mol(NULL), 
	m_SqrMinOverlapFactor(DEFAULT_SQR_MAX_OVERLAP_FACTOR), 
	m_AtomContainerIsVolatile(true),
	m_BondOrder(NULL)
{
	setDefaultPickers();
	setMolecule( _mol );
}

ClashFilterBase::ClashFilterBase( WorkSpace& _mol )
	: FilterBase(), 
	m_Mol(NULL), 
	m_SqrMinOverlapFactor(DEFAULT_SQR_MAX_OVERLAP_FACTOR),
	m_AtomContainerIsVolatile(true),
	m_BondOrder(NULL)
{
	setDefaultPickers();
	setMolecule( _mol );
}

ClashFilterBase::~ClashFilterBase()
{
	// Nothing to clean. CloneHolder deals with all memory management
}

void ClashFilterBase::setOLapFac( double fac )
{
	ASSERT( fac <= 1.0 && fac >= 0.1, ArgumentException, "ClashFilterBase: Call to setOLapFac() fac is invalid value" );
	m_SqrMinOverlapFactor = std::pow(fac,2);
}

void ClashFilterBase::setDefaultPickers()
{
	m_ForPicker = new PickHeavyAtoms(); // CloneHolder deals with memory management
	m_AgainstPicker = new PickHeavyAtoms(); // CloneHolder deals with memory management
}

void ClashFilterBase::setMolecule( MoleculeBase& _mol )
{
	m_Mol = &_mol;
	m_AtomContainerIsVolatile = true; // MoleculeBase atom list can change at any time
	m_BondOrder = NULL; // MolBase has no bondorder
	initCore();
}

void ClashFilterBase::setMolecule( WorkSpace& _mol )
{
	m_Mol = &_mol;
	m_AtomContainerIsVolatile = false; // once allocated, the atoms in a workspace cant change
	m_BondOrder = &_mol.bondorder();
	initCore();
}

void ClashFilterBase::setForPicker( const PickBase& _Picker )
{
	m_ForPicker = _Picker.clone();
	if( m_Mol != NULL ) initCore(); // only init if the molecule has been set, or we will get an exception
}

void ClashFilterBase::setAgainstPicker( const PickBase& _Picker )
{
	m_AgainstPicker = _Picker.clone();
	if( m_Mol != NULL ) initCore(); // only init if the molecule has been set, or we will get an exception
}

void ClashFilterBase::setInternalName()
{
	name = "ClashFilterBase";
}

bool ClashFilterBase::isProxyAtom( int i, int j ) const
{
	if( m_BondOrder != NULL )
	{
		// Use the speedy bondOrder array of the workspace :-D
		return 4 > m_BondOrder->getBondOrder(i,j);
	}
	else
	{
		// We will have to scan the covXXAtom arrays
		ParticleStore& atom = m_Mol->atom;
		for( int q = 0; q < atom[i].cov12atom.size(); q++ )
		{
			if( atom[i].cov12atom[q].i == j ) 
				return true;	
		}
		for( int q = 0; q < atom[i].cov12atom.size(); q++ )
		{
			if( atom[i].cov13atom[q].i == j ) 
				return true;	
		}
		for( int q = 0; q < atom[i].cov12atom.size(); q++ )
		{
			if( atom[i].cov14atom[q].i == j ) 
				return true;	
		}
		return false;
	}
}


ClashFilter::ClashFilter()
	: ClashFilterBase()
{
}

ClashFilter::ClashFilter( MoleculeBase& _mol ) 
	: ClashFilterBase(_mol)
{
}

ClashFilter::ClashFilter( WorkSpace& _mol )
	: ClashFilterBase(_mol)
{
}

ClashFilter::~ClashFilter()
{
	// Nothing, and we like it that way :-D
}

void ClashFilter::initCore()
{
	ASSERT( m_Mol != NULL, CodeException, "Internal NULL reference found in ClashFilterBase. Did you forget to call SetMolecule()?");
	m_ForAtoms.setPicking( *m_Mol, m_ForPicker.data() );
	m_AgainstAtoms.setPicking( *m_Mol, m_AgainstPicker.data() );
}

std::string ClashFilter::reason()
{
	ASSERT( m_Mol != NULL, CodeException, "Internal NULL reference found in ClashFilterBase. Did you forget to call SetMolecule()?");
	
	// Rebuild the cache lists, but only if the container we are using has a volatile atom list.
	if( m_AtomContainerIsVolatile ) initCore(); 

	ParticleStore& atom = m_Mol->atom;
	for( size_t i = 0; i < m_ForAtoms.size(); i++ )
	{
		for( size_t j = 0; j < m_AgainstAtoms.size(); j++ )
		{
			Particle& ai = atom[m_ForAtoms[i]];
			Particle& aj = atom[m_AgainstAtoms[j]];
			if( isProxyAtom(ai.i, aj.i) ) continue; // the same atom will always clash with itself and bonding partners :-D
			double rSqrSum = ai.radius + aj.radius;
			rSqrSum *= rSqrSum;
			double sqrDist = m_Mol->atomxyz(m_ForAtoms[i]).sqrdist(m_Mol->atomxyz(m_AgainstAtoms[j]));
			if( sqrDist < rSqrSum )
			{
				if( (sqrDist / rSqrSum) < m_SqrMinOverlapFactor ) 
				{
					StringBuilder sb;
					sb.setFormat("Clash! Atom(%d,%d) Dist:%8.3fA RadSum:%8.3f OFac:%8.3f")
						(m_ForAtoms[i])(m_AgainstAtoms[j])(sqrt(sqrDist))(sqrt(rSqrSum))(sqrt(m_SqrMinOverlapFactor));
					return sb.toString();
				}
			}
		}
	}
	return std::string("No clashes were detected");
}

bool ClashFilter::passes()
{
	ASSERT( m_Mol != NULL, CodeException, "Internal NULL reference found in ClashFilterBase. Did you forget to call SetMolecule()?");
	
	// The contents of a MoleculeBase can change (are volatile), WorkSpaces are immutable.
	if( m_AtomContainerIsVolatile ) 
		initCore(); // rebuild the cache lists

	ParticleStore& atom = m_Mol->atom;
	for( size_t i = 0; i < m_ForAtoms.size(); i++ )
	{
		Particle& ai = atom[m_ForAtoms[i]];
		for( size_t j = 0; j < m_AgainstAtoms.size(); j++ )
		{			
			Particle& aj = atom[m_AgainstAtoms[j]];
			if( isProxyAtom(ai.i, aj.i) ) continue; // the same atom will always clash with itself and bonding partners :-D
			double rSqrSum = ai.radius + aj.radius;
			rSqrSum *= rSqrSum;
			double sqrDist = m_Mol->atomxyz(m_ForAtoms[i]).sqrdist(m_Mol->atomxyz(m_AgainstAtoms[j]));
			if( sqrDist < rSqrSum )
			{
				if( (sqrDist / rSqrSum) < m_SqrMinOverlapFactor ) 
				{
					return false;
				}
			}
		}
	}
	return true;
}

const double DEFAULT_GRID_CUTOFF = 5.0;

ClashGridFilter::ClashGridFilter()
	: ClashFilterBase()
{
	GridCutoff = DEFAULT_GRID_CUTOFF;
}

ClashGridFilter::ClashGridFilter( MoleculeBase& _mol ) 
	: ClashFilterBase(_mol)
{
	GridCutoff = DEFAULT_GRID_CUTOFF;
}

ClashGridFilter::ClashGridFilter( WorkSpace& _mol )
	: ClashFilterBase(_mol)
{
	GridCutoff = DEFAULT_GRID_CUTOFF;
}

ClashGridFilter::~ClashGridFilter()
{
	// Nothing, and we like it that way :-D
}

void ClashGridFilter::initCore()
{
	ASSERT( m_Mol != NULL, CodeException, "Internal NULL reference found in ClashGridFilter. Did you forget to call SetMolecule()?");
	m_ForAtoms.setPicking( *m_Mol, m_ForPicker.data() );
	m_AgainstAtoms.setPicking( *m_Mol, m_AgainstPicker.data() );
	m_ProxyGrid = ProximityGrid( *m_Mol, m_AgainstPicker.data(), GridCutoff );
}
	
bool ClashGridFilter::passes()
{
	ASSERT( m_Mol != NULL, CodeException, "Internal NULL reference found in ClashGridFilter. Did you forget to call SetMolecule()?");

	// The contents of a MoleculeBase can change (are volatile), WorkSpaces are immutable.
	if( m_AtomContainerIsVolatile ) 
		initCore(); // rebuild the cache lists

	// Internally calls refresh()
	m_ProxyGrid->setCutoff( GridCutoff );

	ASSERT( &m_ProxyGrid.data() != NULL, CodeException, "Proxy grid has not been created!!");

	const GridPoint* gridPoint = NULL;
	ParticleStore& atom = m_Mol->atom;
	for( size_t i = 0; i < m_ForAtoms.size(); i++ )
	{
		Particle& ai = atom[m_ForAtoms[i]];
		if( m_ProxyGrid->atomList( ai.pos(), gridPoint ) )
		{
			for( size_t j = 0; j < gridPoint->size(); j++ )
			{				
				Particle& aj = atom[gridPoint->atomIndexes[j]];
				if( isProxyAtom(ai.i, aj.i) ) 
					continue; // the same atom will always clash with itself and bonding partners :-D
				double rSqrSum = ai.radius + aj.radius;
				rSqrSum *= rSqrSum;
				double sqrDist = m_Mol->atomxyz(ai.i).sqrdist(m_Mol->atomxyz(aj.i));
				if( sqrDist < rSqrSum )
				{
					if( (sqrDist / rSqrSum) < m_SqrMinOverlapFactor ) 
					{
						return false;
					}
				}
			}
		}
	}
	return true;
}

std::string ClashGridFilter::reason()
{
	ASSERT( m_Mol != NULL, CodeException, "Internal NULL reference found in ClashGridFilter. Did you forget to call SetMolecule()?");

	// The contents of a MoleculeBase can change (are volatile), WorkSpaces are immutable.
	if( m_AtomContainerIsVolatile ) 
		initCore(); // rebuild the cache lists

	m_ProxyGrid->setCutoff( GridCutoff );

	ASSERT( &m_ProxyGrid.data() != NULL, CodeException, "Proxy grid has not been created!!");

	const GridPoint* gridPoint = NULL;
	ParticleStore& atom = m_Mol->atom;
	for( size_t i = 0; i < m_ForAtoms.size(); i++ )
	{
		Particle& ai = atom[m_ForAtoms[i]];
		if( m_ProxyGrid->atomList( ai.pos(), gridPoint ) )
		{
			for( size_t j = 0; j < gridPoint->size(); j++ )
			{				
				Particle& aj = atom[gridPoint->atomIndexes[j]];
				if( isProxyAtom(ai.i, aj.i) ) 
					continue; // the same atom will always clash with itself and bonding partners :-D
				double rSqrSum = ai.radius + aj.radius;
				rSqrSum *= rSqrSum;
				double sqrDist = m_Mol->atomxyz(ai.i).sqrdist(m_Mol->atomxyz(aj.i));
				if( sqrDist < rSqrSum )
				{
					if( (sqrDist / rSqrSum) < m_SqrMinOverlapFactor ) 
					{
						StringBuilder sb;
						sb.setFormat("Clash! Atom(%d,%d) Dist:%8.3fA RadSum:%8.3f OFac:%8.3f")
							(ai.i)(aj.i)(sqrt(sqrDist))(sqrt(rSqrSum))(sqrt(m_SqrMinOverlapFactor));
						return sb.toString();
					}
				}
			}
		}
	}
	return std::string("No clashes were detected");
}

