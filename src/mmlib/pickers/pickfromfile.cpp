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
#include "system/fundamentals.h"
#include "system/molecule.h"
#include "pickfromfile.h"

ResidueList::ResidueList()
{
}

void ResidueList::setChainID(const char& _chainID)
{
	m_ChainID = _chainID;
}

char ResidueList::getChainID() const
{
	return m_ChainID;
}

void ResidueList::pushback(const Sequence::BioSource& _resInfo)
{
	m_ResList.push_back(_resInfo);
}

int ResidueList::getResListSize() const
{
	return m_ResList.size();
}

int ResidueList::getBiologicalIndex(int _ResListIndex) const
{
	return m_ResList[_ResListIndex].getBiologicalIndex();
}

char ResidueList::getBiologicalIcode(int _ResListIndex) const
{
	return m_ResList[_ResListIndex].getBiologicalICode();
}

std::string ResidueList::getResName(int _ResListIndex) const
{
	return m_ResList[_ResListIndex].getResName();
}

PickResiduesFromFile::PickResiduesFromFile(const MoleculeBase& _molBase)
	:PickResidueBase()
{
	// Nothing happens when this default constructor is called.
	// matches() will therefore return false for all particles/residues.
}

PickResiduesFromFile::PickResiduesFromFile(const MoleculeBase& _molBase, const std::string& _filename)
	:PickResidueBase()
{
	loadList(_filename);
	compareAll(_molBase);
}

PickResiduesFromFile* PickResiduesFromFile::clone() const
{
	return new PickResiduesFromFile(*this);
}

void PickResiduesFromFile::loadList(const std::string& _filename)
{
	std::string _chainIDString = " "; ///< local store of chainID as a string.
	char _chainID = CHAR_MAX;         ///< local store of chainID as a char.
	std::string _iresString = " ";    ///< local store of the biological residue index as a string.
	int _ires = -1;                   ///< local store of the biological residue index as an int.
	std::string _residueName = " ";   ///< local store of residue name (3 letter code).

	const char *whiteSpace = " \t\12\15"; ///< Whitespace can be "space", tab, newline character or end-of-line character
	const char *commentChar = "#\12\15";  ///< Possible comment characters

	std::ifstream rotchoiceFile(_filename.c_str(), std::ifstream::in);
	if( !rotchoiceFile.is_open() ) throw IOException("Rotamer choices file " + _filename + " not found." );

	std::string line;

	while(std::getline(rotchoiceFile, line))
	{
		// _IcodeString needs to be initialised here, within the while loop, so that it defaults to a space.
		std::string _IcodeString = " ";   ///< local store of Icode as a string. If there's no Icode specified then store a space.
		char _Icode = ' ';                ///< local store of Icode as a char. If there's no Icode specified then store a space.

		// Reset local variables to allow error checking.
		_chainIDString = " ";
		_chainID = CHAR_MAX;
		_iresString = " ";
		_ires = -1;
		_residueName = " ";

		removecomments(line, commentChar);
		std::vector<std::string> token;

		token = chopstr(line, whiteSpace); // each whitespace-separated value on the line will become a member of the array.
		if(token.size() <= 0) continue; // this will be an empty line
		if(token.size() > 4)            // each line should only contain a chain ID, a residue number, a residue 3 letter code and a residue Icode (optional).
		{
			throw ProcedureException("Each line in Rotamer choice file should only include:\n\ta chain ID, a residue number, a residue 3 letter code (or *) and a residue Icode (optional).");
		}

		_chainIDString = token[0]; // since a string is an array of chars the first char in the array is now the chainID.
		_chainID = _chainIDString[0];
		_iresString = token[1];    // don't store this...store an int, _ires, instead
		str2int(_iresString,_ires);
		_residueName = token[2];

		if(token.size() == 4)
		{
			_IcodeString = token[3]; // since a string is an array of chars the first char in the array is now the Icode.
			_Icode = _IcodeString[0];
		}

		Sequence::BioSource _tempResInfo(_ires,_Icode,_residueName); ///< temporary local store of residue info before adding it to m_ResListByChain.

		// Error checking - if these variables are still at their starting values throw an exception to warn the user.
		if( _chainID != CHAR_MAX
			&& _tempResInfo.getBiologicalIndex() != - 1
			&& !(cmpstring(_tempResInfo.getResName()," ")) )
		{
			// If any residues have already been added to the list (m_ResListByChain) then go through the list and
			// check whether there is already a ResidueList with the chainID of the current residue. If there is,
			// add the current residue to that ResidueList rather than creating a new one for it. This way, residues
			// will be grouped by chainID.
			if( m_ResListByChain.size() > 0 )
			{
				// loop over the list
				for( int i = 0; i < m_ResListByChain.size(); i++ )
				{
					if( _chainID == m_ResListByChain[i].getChainID() )
					{
						// if that chainID already exists then just add the BioSource information for this residue
						m_ResListByChain[i].pushback(_tempResInfo);
						break;
					}

					// this else should be true only if the for loop has finished and the previous if statement hasn't been hit.
					else if( i == ( m_ResListByChain.size() -1 ) )
					{
						// if that chainID doesn't already exist then make a whole new ResidueList containing the
						// chainID and BioSource information for the current residue.
						ResidueList _tempResList; ///< temporary local store of residue info and chainID before adding it to m_ResListByChain.
						_tempResList.setChainID(_chainID);
						_tempResList.pushback(_tempResInfo);
						m_ResListByChain.push_back(_tempResList);
						break;
					}
				}
			}

			// If nothing has been added to m_ResListByChain yet then make a whole new ResidueList containing the
			// chainID and BioSource information for the current residue.
			else if( m_ResListByChain.size() == 0 )
			{
				ResidueList _tempResList; ///< temporary local store of residue info and chainID before adding it to m_ResListByChain.
				_tempResList.setChainID(_chainID);
				_tempResList.pushback(_tempResInfo);
				m_ResListByChain.push_back(_tempResList);
			}
		}
		else
			throw ParseException("Rotamer choice file set incorrectly:\neach line (residue choice) must include a chain ID, a residue number, a residue 3 letter code (or *) and a residue Icode (optional).");
	}
}

bool PickResiduesFromFile::matches( const Particle& particle ) const
{
	size_t _iresSequential = particle.ir;
	return m_ResIsPicked[_iresSequential];
}

bool PickResiduesFromFile::matches( const Residue& _res ) const
{
	size_t _iresSequential = _res.ir;
	return m_ResIsPicked[_iresSequential];
}

void PickResiduesFromFile::compareAll(const MoleculeBase& _molBase)
{
	size_t _molBaseNRes = _molBase.nResidues(); ///< number of residues in the "MoleculeBase" (probably WorkSpace).

	m_ResIsPicked.clear();
	m_ResIsPicked.resize(_molBaseNRes,false); ///< Size is set to the number of residues in the "_molecule" and each bool is set to false.

	int _sizeOfResListByChain = m_ResListByChain.size(); // local store of the number of chains

	// loop over each chain (picked residues are grouped by chain)
	for(int i = 0; i < _sizeOfResListByChain; i++)
	{
		char _chainID_i = m_ResListByChain[i].getChainID(); // what chain ID does this group have?

		// loop over all residues in the "MoleculeBase" (probably WorkSpace).
		for(int ir = 0; ir < _molBaseNRes; ir++)
		{
			if( m_ResIsPicked[ir] ) continue; ///< If the residue has already been picked then continue to the next loop.

			Residue _molBaseRes = _molBase.getRes(ir);                     // local store of the current residue
			int _molBaseResIFirst = _molBaseRes.ifirst;                    // the first atom of this residue
			char _molBaseChainID = _molBase.getChainID(_molBaseResIFirst); // the chainID of that atom (since Residue doesn't know it's chainID).

			// if this residue doesn't have an ifirst
			if( _molBaseResIFirst == -1 )
			{
				StringBuilder sb;
				sb.setFormat("Residue %d does not define iFirst")(ir);
				THROW(ProcedureException,sb.toString());
			}

			if( _chainID_i != _molBaseChainID ) continue; ///< if chainID isn't the same then it's not the same residue.

			// loop over all residues picked for this chain
			for(int j = 0; j < m_ResListByChain[i].getResListSize(); j++)
			{
				//Sequence::ResidueInfo _checkResInfo = _molBase.getSequence().getResidue(ir);
				if( m_ResListByChain[i][j].isEqual(_molBase.getSequence().getResidue(ir)) )
				{
					// if this residue (in m_ResListByChain[i].m_ResList[j]) is the same as the current residue in
					//the "MoleculeBase" (probably WorkSpace) then set that "MoleculeBase" residue as picked.
					m_ResIsPicked[ir] = true;
				}
			}
		}
	}
}




