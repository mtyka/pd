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
#include <fstream>
#include "tools/stringtool.h"
#include "pickers/basicpickers.h"
#include "workspace/workspace.h"
#include "workspace/bondorder.h"
#include "pops.h"

void PopsAtomDBType::parse( const std::string& _Line )
{
	std::vector<std::string> parts = chopstr( _Line, TOKEN_WHITESPACE.c_str() );
	if( parts.size() != 6 ) throw ParseException("POPS: Atom set initial line contains the wrong number of elements");
	if( 0 != parts[0].compare("ATOM") ) THROW( CodeException, "Only atom lines should ever be fed in here!" );

	atomName = parts[1];
	resName = parts[2];

	int endCol;
	if( 
		0 != str2double( parts[3], radius ) ||
		0 != str2double( parts[4], param ) ||
		0 != str2int( parts[5], endCol )
		) throw ParseException("POPS: Invalid value in ATOM line");

	hydrophilic = (endCol != 0);

	return;
}

void PopsDat::parse( const std::string& _Line )
{
	std::vector<std::string> parts = chopstr( _Line, TOKEN_WHITESPACE.c_str() );
	if( parts.size() != 5 ) throw ParseException("POPS: Atom set initial line contains the wrong number of elements");
	if( 
		0 != str2int( parts[0], expectedAtomTypeCount ) ||
		0 != str2double( parts[1], b12 ) ||
		0 != str2double( parts[2], b13 ) ||
		0 != str2double( parts[3], b14 ) ||
		0 != str2double( parts[4], bOther )
		) throw ParseException("POPS: Invalid value in atom set initial line");
	return;
}

Pops::Pops() 
	: m_WSpace(NULL), 
	m_CurrentDat(NULL),
	m_Mode(AllAtom),
	ProbeRadius(1.4),
	m_DataRead(false)
{
}

// Stand-alone helper parsing function
inline void obtainLine( std::ifstream& _stream, std::string& line )
{
	// Pump lines out of the file, discarding those which are comments and empty
	// Verify we are not at EOF at each case
	while(std::getline(_stream,line))
	{		
		if( _stream.fail() ) throw ParseException("POPS: Unexpected 'stream error' whilst parsing data");	
		line = trim( line, TOKEN_WHITESPACE );
		if( line.size() == 0 ||  line[0] == '#' )
		{
			if( _stream.eof() ) throw ParseException("POPS: Unexpected 'end of file' whilst parsing data");
			continue; // ignore blank lines and comments
		}			
		break;
	}
}

void Pops::readDat( const std::string& _DatFileName )
{
	std::ifstream file(_DatFileName.c_str(), std::ifstream::in);

	if( !file.is_open() ) throw(IOException( "POPS: Data file not found: '" + _DatFileName + "'!" ));

	std::string line; // Temporary line container

	// Find file data begining, assert that this is a POPS dat file
	while(true)
	{
		obtainLine(file,line);
		if( 0 == line.compare("BEGIN POPS") )
		{
			break;
		}
		throw ParseException("POPS: POPS dat format is not correct");
	}

	PopsDat* dat = NULL;
	// Extract our data
	while(true)
	{
		obtainLine(file,line);
		if( 0 == line.compare(0,4,"ATOM") )
		{
			if( dat == NULL )
			{
				throw ParseException("POPS: Atoms are being defined without a corresponding BEGIN statment.");
			}
			PopsAtomDBType a;
			a.parse( line );
			dat->atoms.push_back(a);
		}
		else if( 0 == line.compare("BEGIN ATOMIC") )
		{
			// We are defining an atomic resolution set
			// clear all existing data
			m_Atomic = PopsDat();
			dat = &m_Atomic;
			obtainLine(file,line);
			dat->parse(line);
		}
		else if( 0 == line.compare("BEGIN COARSE") )
		{
			// We are defining a coarsse resolution set
			// clear all existing data
			m_Coarse = PopsDat();
			dat = &m_Coarse;
			obtainLine(file,line);
			dat->parse(line);
		}
		else if( 0 == line.compare("END") )
		{
			if( dat == NULL ) 
			{
				// This is the file's concluding END statement
				break;
			}
			else
			{
				// We are ending a BEGIN xxx statement
				if( dat->expectedAtomTypeCount != (int)dat->atoms.size() )
				{
					throw ParseException("POPS: The incorrect number of atom have been imported.");
				}
				dat = NULL;
			}
		}
		else
		{
			printf("POPS: Parse WARNING: Uninterpreted line '%s'\n", line.c_str() );
		}
	}
	m_DataRead = true;

	file.close();
}

void Pops::setTo( const WorkSpace& _WSpace, PopsMode _Mode )
{
	m_WSpace = &_WSpace;
	m_Mode = _Mode;

	m_AtomIndexes.clear();
    m_Atoms.clear();

	switch( m_Mode )
	{
	case AllAtom:
		{
			m_CurrentDat = &m_Atomic;
			if( m_Atomic.atoms.size() == 0 )
			{
				throw ProcedureException("POPS: No data has yet been loaded for AllAtom mode.");
			}
			for( size_t i = 0; i < _WSpace.nAtoms(); i++ )
			{
				if( _WSpace.atom[i].isHydrogen() )
				{
					m_AtomIndexes.push_back(SIZE_T_FAIL);
					continue;
				}
				else
				{
					m_AtomIndexes.push_back(m_Atoms.size());
				}

				bool found = false;
				for( size_t j = 0; j < m_Atomic.atoms.size(); j++ )
				{
					if( 0 == m_Atomic.atoms[j].atomName.compare( _WSpace.atom[i].pdbname ) &&
						0 == m_Atomic.atoms[j].resName.compare( _WSpace.atom[i].parentl3name ) )
					{
						m_Atoms.push_back( PopsAtom( m_Atomic.atoms[j], i, &_WSpace.atomxyz(i)) );
						found = true;
						break;
					}
				}
				if( !found )
				{
					for( size_t j = 0; j < m_Atomic.atoms.size(); j++ )
					{
						if( 0 == m_Atomic.atoms[j].resName.compare( "*" ) &&
							0 == m_Atomic.atoms[j].atomName.compare( _WSpace.atom[i].pdbname ) )
						{
							m_Atoms.push_back( PopsAtom( m_Atomic.atoms[j], i, &_WSpace.atomxyz(i)) );
							found = true;
							break;
						}
					}
				}
				if( !found )
				{
					StringBuilder sb;
					sb.setFormat("POPS: Cannot find atom type %s %s in the database.\n")
						(_WSpace.atom[i].pdbname)
						(_WSpace.atom[i].parentl3name);
					throw ProcedureException(sb.toString());
				}
			}
			break;
		}
	case Coarse:
		{
			m_CurrentDat = &m_Coarse;
			if( m_Coarse.atoms.size() == 0 )
			{
				throw ProcedureException("POPS: No data has yet been loaded for AllAtom mode.");
			}
			for( size_t i = 0; i < _WSpace.nAtoms(); i++ )
			{
				if( !_WSpace.atom[i].isCAlpha() )
				{
					m_AtomIndexes.push_back(SIZE_T_FAIL);
					continue;
				}
				else
				{
					m_AtomIndexes.push_back(m_Atoms.size());
				}

				bool found = false;
				for( size_t j = 0; j < m_Coarse.atoms.size(); j++ )
				{
					if( 0 == m_Coarse.atoms[j].atomName.compare( _WSpace.atom[i].pdbname ) &&
						0 == m_Coarse.atoms[j].resName.compare( _WSpace.atom[i].parentl3name ) )
					{
						m_Atoms.push_back( PopsAtom( m_Coarse.atoms[j], i, &_WSpace.atomxyz(i)) );
						found = true;
						break;
					}
				}
				if( !found )
				{
					for( size_t j = 0; j < m_Coarse.atoms.size(); j++ )
					{
						if( 0 == m_Coarse.atoms[j].resName.compare( "*" ) &&
							0 == m_Coarse.atoms[j].atomName.compare( _WSpace.atom[i].pdbname ) )
						{
							m_Atoms.push_back( PopsAtom( m_Coarse.atoms[j], i, &_WSpace.atomxyz(i) ) );
							found = true;
							break;
						}
					}
				}
				if( !found )
				{
					StringBuilder sb;
					sb.setFormat("POPS: Cannot find atom type %s %s in the database.\n")
						(_WSpace.atom[i].pdbname)
						(_WSpace.atom[i].parentl3name);
					throw ProcedureException(sb.toString());
				}				
			}
			break;
		}
	default:
		{
			THROW( CodeException, "Unknown PopsMode encountered!");
		}
	}
}

void Pops::CoreAsserts() const
{
	ASSERT( m_DataRead, CodeException, "POPS::readDat() has not been called" );
	ASSERT( m_WSpace != NULL, CodeException, "POPS::setTo() has not been called" );
	D_ASSERT( m_WSpace->nAtoms() == m_AtomIndexes.size(), CodeException, "Internal Code Error");
}

void Pops::info() const
{
	CoreAsserts();
	OverallSASA sasa = getSASA();

	switch( m_Mode )
	{
	case AllAtom:
		{
			printf("Pops All-Atom Info:\n");
			break;
		}
	case Coarse:
		{
			printf("Pops Coarse Info:\n");
			break;
		}
	default:
		{
			THROW( CodeException, "Unknown PopsMode encountered!");
		}
	}	

	printf("\tSASA atoms: %d\n", m_Atoms.size() );
	printf("\tProbe Radius: %5.2lf\n", ProbeRadius );
	printf("\tTotal SASA: %8.3lf\n", sasa.SASA );
	printf("\tHydrophobic SASA: %8.3lf\n", sasa.hydrophobicSASA );
	printf("\tHydrophilic SASA: %8.3lf\n", sasa.hydrophilicSASA );
}

void Pops::detail() const
{
	CoreAsserts();
	info();
	const ParticleStore& atom = m_WSpace->atom;
	printf("\nSASA Particle List:\n");
	for( size_t i = 0; i < m_Atoms.size(); i++ )
	{
		atom[m_Atoms[i].parentIndex].info(Verbosity::Normal,4,false);
		printf(" Radius: %6.3lf SASA: %6.3lf NOverlap %d\n", m_Atoms[i].radius, m_Atoms[i].sasa, m_Atoms[i].NOverlap );
	}
}

OverallSASA Pops::getSASA() const
{
	CoreAsserts();
	return m_SASAInfo;
}

double Pops::sasaFraction( const PickAtomRange& _Range ) const
{
	size_t startIndex = _Range.getStartAtomIndex();
	size_t endIndex = _Range.getEndAtomIndex();
	double sasa = 0.0;
	double maxSASA = 0.0;
	for( size_t i = startIndex; i <= endIndex; i++ )
	{
		size_t index = m_AtomIndexes[i];
		if( index != SIZE_T_FAIL )
		{
			// This is in the include list dependent on mode and atom selection
			sasa += m_Atoms[index].sasa;
			maxSASA += m_Atoms[index].maxSASA;
		}
	}
	return sasa / maxSASA;
}

double Pops::SASA( const PickAtomRange& _Range ) const
{
	size_t startIndex = _Range.getStartAtomIndex();
	size_t endIndex = _Range.getEndAtomIndex();
	double sasa = 0.0;
	for( size_t i = startIndex; i <= endIndex; i++ )
	{
		size_t index = m_AtomIndexes[i];
		if( index != SIZE_T_FAIL )
		{
			// This is in the include list dependent on mode and atom selection
			sasa += m_Atoms[index].sasa;
		}
	}
	return sasa;
}

double Pops::atomSASA( int ia ) const
{
	CoreAsserts();
	ASSERT( ia >=0 && ia < m_WSpace->nAtoms(), OutOfRangeException, "Atom request is out of range");
	size_t index = m_AtomIndexes[ia];
	return index == SIZE_T_FAIL ? 0.0 : m_Atoms[index].sasa;
}
double Pops::resFraction( int ir ) const
{
	CoreAsserts();
	ASSERT( ir >=0 && ir < m_WSpace->nResidues(), OutOfRangeException, "Atom request is out of range");
	switch( m_Mode )
	{
	case AllAtom:
		{
			int resStart =  m_WSpace->res[ir].ifirst;
			int resEnd =  m_WSpace->res[ir].ilast;
			ASSERT( resStart != -1 && resEnd != -1, OutOfRangeException, "Residue range is not defined");
			double sum = 0.0f;
			double sumMax = 0.0f;
			for( int i = resStart; i <= resEnd; i++ )
			{
				size_t lookup = m_AtomIndexes[i];
				if( lookup != SIZE_T_FAIL )
				{
					double sasa = m_Atoms[lookup].sasa;
					D_ASSERT( sasa != -1, CodeException, "Assumption is not true");
					sum += sasa;
					sumMax += m_Atoms[lookup].maxSASA;
				}
			}
			return sum / sumMax;
		}
	case Coarse:
		{
			int CAIndex = m_WSpace->res[ir].iCA;
			ASSERT( CAIndex != -1, ProcedureException, "CA Position is undefined in the residue");
			size_t index = m_AtomIndexes[CAIndex];
			ASSERT( m_AtomIndexes[CAIndex] != SIZE_T_FAIL, CodeException, "Assumption is not true");			
			return m_Atoms[index].sasa / m_Atoms[index].maxSASA;
		}
	default:
		{
			THROW( CodeException, "Unknown PopsMode encountered!");
		}
	}
}

double Pops::resSASA( int ir ) const
{
	CoreAsserts();
	ASSERT( ir >=0 && ir < m_WSpace->nResidues(), OutOfRangeException, "Atom request is out of range");
	switch( m_Mode )
	{
	case AllAtom:
		{
			int resStart =  m_WSpace->res[ir].ifirst;
			int resEnd =  m_WSpace->res[ir].ilast;
			ASSERT( resStart != -1 && resEnd != -1, OutOfRangeException, "Residue range is not defined");
			double sum = 0.0f;
			for( int i = resStart; i <= resEnd; i++ )
			{
				size_t lookup = m_AtomIndexes[i];
				if( lookup != SIZE_T_FAIL )
				{
					double sasa = m_Atoms[lookup].sasa;
					D_ASSERT( sasa != -1, CodeException, "Assumption is not true");
					sum += sasa;
				}
			}
			return sum;
		}
	case Coarse:
		{
			int CAIndex = m_WSpace->res[ir].iCA;
			ASSERT( CAIndex != -1, ProcedureException, "CA Position is undefined in the residue");
			size_t index = m_AtomIndexes[CAIndex];
			ASSERT( m_AtomIndexes[CAIndex] != SIZE_T_FAIL, CodeException, "Assumption is not true");			
			return m_Atoms[index].sasa;
		}
	default:
		{
			THROW( CodeException, "Unknown PopsMode encountered!");
		}
	}
}

void Pops::calc()
{
	CoreAsserts();

	ASSERT( ProbeRadius > 1.0 && ProbeRadius < 2.0, ArgumentException, "Pops: 'ProbeRadius' (H20) has been set to a non-sensical value.");

	// Proxy
	const ParticleStore& atom = m_WSpace->atom;
	const BondOrder& bond = m_WSpace->bondorder();
	size_t nSASAAtom = m_Atoms.size();

	m_SASAInfo = OverallSASA(); // Reset to 0.0's
	double twoProbe = 2.0 * ProbeRadius; // Initialise 'double the probe radius' - ProbeRadius may have been changed.

	// Calculate all atomic SASAs
	for( size_t i = 0; i < nSASAAtom; i++ )
	{		
		// Atom i
		PopsAtom& si = m_Atoms[i];
		int qi = si.parentIndex;
		double Ri = si.radius;

		// calc the 'Si' and its inverse for this atom
		double Si = si.radius + ProbeRadius;
		Si *= Si; // square it
		Si *= Maths::MathConst::FourPI;
		double invSi = 1.0 / Si; // the inverse is used repeatedly below

		// Now lets calc the sasa
		si.NOverlap = 0;	
		si.sasa = Si; // Initialise to the max possible SASA
		si.maxSASA = Si; // This is the maxSASA possible

		// Loop over other atoms in the list
		for( size_t j = 0; j < nSASAAtom; j++ )
		{
			if( i == j ) continue; // You cant occlude SASA from yourself!

			const PopsAtom& sj = m_Atoms[j];
			int qj = sj.parentIndex;
			double Rj = sj.radius;
			
			double rij = si.pos->dist( *sj.pos );
			double cutoff = Rj + Ri + twoProbe;		
			if( rij < cutoff )
			{
				si.NOverlap++; // We have overlapping spheres

				int orderij;
				if( m_Mode == Coarse )
				{
					// When we are talking about just CA atoms:
					// order takes on new meaning as the difference in residue index
					orderij = abs(atom[qi].ir-atom[qj].ir);	
					// BTW - The official POPS cannot do this, and seems to detect bond order by using 'rij' -
					// doing something like this: 'if( rij > Ri ) orderij = 4;'.
					// I found this by putting loads of test cases into their webserver. 
					// My test case was to take two ALA CA atoms and gradually move them appart and measure SASA
					// using POPS and PD. The POPS algorith yields a large step in SASA when rij gets above Ri, as the
					// detection changes the bond order to a 1-4 pair.
					// Basically for odd cases, their bondorder is 'wrong', because of this detection. 
					// For protein normal examples, both implementations match exactly.
				}
				else
				{
					orderij = bond.getBondOrder( qi, qj );				
				}

				// Get the paramter pij dependent on the bond order
				double pij;
				if( orderij >= 4 )
				{
					pij = m_CurrentDat->bOther;
				}
				else if( orderij == 3 )
				{
					pij = m_CurrentDat->b14;
				}
				else if( orderij == 2 )
				{
					pij = m_CurrentDat->b13;
				}
				else // orderij == 1, cant be 0 
				{
					D_ASSERT( orderij != 0, CodeException, "Pops: Internal bond order error, bond order cannot be '0' if i != j");
					pij = m_CurrentDat->b12;
				}

				double bij = Maths::MathConst::PI * 
					(Ri + ProbeRadius) *
					(Ri + Rj + twoProbe - rij) * 
					(1.0+((Rj-Ri)*(1.0/rij)));
				si.sasa *= 1.0 - (si.param * pij * bij * invSi);
			}
		}

		if( m_Atoms[i].hydrophilic )
		{
			m_SASAInfo.hydrophilicSASA += m_Atoms[i].sasa;
		}
		else
		{
			m_SASAInfo.hydrophobicSASA += m_Atoms[i].sasa;
		}

		m_SASAInfo.SASA += m_Atoms[i].sasa;
	}
}




