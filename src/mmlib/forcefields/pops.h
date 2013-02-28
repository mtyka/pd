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

#ifndef __POPS_H
#define __POPS_H

//   1. Fraternali, F. and Cavallo, L.
//      Parameter optimized surfaces (POPS): analysis of key interactions and conformational changes in the ribosome.
//      Nucleic Acids Research 30 (2002) 2950-2960.
//
//   2. Cavallo, L., Kleinjung, J. and Fraternali, F.
//      POPS: A fast algorithm for solvent accessible surface areas at atomic and residue level.
//      Nucleic Acids Res. 31 (2003) 3364-3366.
//
//   3. Kleinjung, J. and Fraternali, F.
//      POPSCOMP: an automated interaction analysis of biomolecular complexes.
//      Nucleic Acids Research 33 (2005) W342-W346.

#include <string>
#include "workspace/workspace.fwd.h"

// Forward declarations
class PickAtomRange;

struct PopsAtomBase
{
	PopsAtomBase() : radius(0.0), param(0.0), hydrophilic(false)
	{
	}
	double radius;
	double param;
	bool hydrophilic;
};

struct PopsAtomDBType : public PopsAtomBase
{
	PopsAtomDBType()
	{
	}

	void parse( const std::string& _Line );

	std::string atomName;
	std::string resName;
};

struct PopsAtom : public PopsAtomBase
{
	PopsAtom() : sasa(0.0), parentIndex(-1), NOverlap(0), pos(NULL), maxSASA(0.0)
	{
	}

	PopsAtom( const PopsAtomDBType& cloneLib, int _ParentIndex, const Maths::dvector* _pos )
		: sasa(-1.0),
		NOverlap(0),
		pos( _pos ),
		maxSASA(0.0)
	{
		radius = cloneLib.radius;
		param = cloneLib.param;
		hydrophilic = cloneLib.hydrophilic;
		parentIndex = _ParentIndex;
	}

	int parentIndex;
	int NOverlap;
	double sasa;
	double maxSASA;
	const Maths::dvector* pos;
};

struct PopsDat
{
	PopsDat() : expectedAtomTypeCount(0), b12(0.0), b13(0.0), b14(0.0), bOther(0.0)
	{
	}

	void parse( const std::string& _Line );

	int expectedAtomTypeCount;
	std::vector<PopsAtomDBType> atoms;
	double b12;
	double b13;
	double b14;
	double bOther;
};

struct OverallSASA
{
	OverallSASA()
	{
		SASA = 0.0;
		hydrophilicSASA = 0.0;
		hydrophobicSASA = 0.0;
	}

	double SASA;
	double hydrophilicSASA;
	double hydrophobicSASA;
};






//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///   1. Fraternali, F. and Cavallo, L.
///      Parameter optimized surfaces (POPS): analysis of key interactions and conformational changes in the ribosome.
///      Nucleic Acids Research 30 (2002) 2950-2960.
///
///   2. Cavallo, L., Kleinjung, J. and Fraternali, F.
///      POPS: A fast algorithm for solvent accessible surface areas at atomic and residue level.
///      Nucleic Acids Res. 31 (2003) 3364-3366.
///
///   3. Kleinjung, J. and Fraternali, F.
///      POPSCOMP: an automated interaction analysis of biomolecular complexes.
///      Nucleic Acids Research 33 (2005) W342-W346.
///
/// \author  Jon Rea 
///
/// \todo 
///
/// \bug 
///
class PD_API Pops
{
public:
	enum PopsMode
	{
		Coarse,
		AllAtom
	};

	Pops();

	// Setup
	void readDat(const std::string& _DatFileName);
	void setTo( const WorkSpace& _WSpace, PopsMode _Mode = AllAtom );

	// Info
	void detail() const;
	void info() const;

	/// Calculate the SASA
	void calc();

	// Query
	OverallSASA getSASA() const;
	double sasaFraction( const PickAtomRange& _Range ) const;
	double SASA( const PickAtomRange& _Range ) const;
	double atomSASA( int ia ) const;
	double resSASA( int ir ) const;
	double resFraction( int ir ) const;

	// Public settings
	double ProbeRadius;

private:
	// Private functions
	inline void CoreAsserts() const;

	// Internal flags
	bool m_DataRead;
	PopsMode m_Mode;

	// Member data and proxies
	const WorkSpace* m_WSpace;
	OverallSASA m_SASAInfo;
	PopsDat* m_CurrentDat;
	PopsDat m_Coarse;
	PopsDat m_Atomic;
	std::vector<size_t> m_AtomIndexes;
	std::vector<PopsAtom> m_Atoms;
};

#endif


