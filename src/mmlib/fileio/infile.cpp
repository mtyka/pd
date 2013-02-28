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

#include <iostream>

#include "sequence/sequence.h"
#include "sequence/alignment.h"
#include "system/molecule.h"

#include "infile.h"

namespace IO
{
	void PD_API FileParticle::info() const
	{
		std::cout << ' ' <<  atomNum
			<< ' ' <<  atomName
			<< ' ' <<  resName
			<< ' ' <<  chainID
			<< ' ' <<  iCode
			<< ' ' <<  resNum
			<< ' ' <<  x
			<< ' ' <<  y
			<< ' ' <<  z
			<< std::endl;
	}

	void setPolymerPositions( const FFParamSet& ffps, MoleculeBase& mol, const Sequence::AlignmentDef& sequencePair, const std::vector<FileParticle>& particles, Verbosity::Type verbose )
	{
		ASSERT(sequencePair.isSeqAligned(),ProcedureException,"ExpSeqPair is not aligned!");

		const int* equiv = sequencePair.getEquivelency();
		const Sequence::BioSequence& bioseq = sequencePair.getSeq1();
		const Sequence::BioSequence& strucseq = sequencePair.getSeq2();
		const Sequence::BioSequence& molseq = mol.getSequence();
		ParticleStore& atom = mol.atom;

		// The molecule *SHOULD* have been made from the biological sequence if this function is being called
		ASSERT(bioseq==molseq,ProcedureException,"MoleculeBase:ExpSeqPair mapping appears incorrect!");

		// check we have something to work with...
		if( particles.size() == 0 || mol.nResidues() == 0 ) return;

		// Detect residue boundries in the particle array
		int resSeqIndex = 0; // The 'particles' have no internal residue serial indexer, we have to count it here...
		char previousICode = particles[0].iCode;
		std::string previousName = particles[0].resName;
		int previousResNum = particles[0].resNum;
		std::vector<int> particleStart;
		std::vector<int> particleEnd;
		particleStart.push_back(0);
		for(int i = 0; i < particles.size(); i++)
		{
			// find atoms in the user-supplied PDB file
			if( (particles[i].iCode != previousICode) ||
				(particles[i].resNum != previousResNum)||
				(particles[i].resName.compare(previousName) != 0 ) )
			{
				particleEnd.push_back(i-1);
				particleStart.push_back(i);
				previousICode = particles[i].iCode;
				previousResNum = particles[i].resNum;
				previousName = particles[i].resName;
				resSeqIndex++;
			}
		}
		resSeqIndex++;
		particleEnd.push_back((int)particles.size()-1);
		ASSERT(strucseq.size() == resSeqIndex, CodeException, "Assumption failed");	

		// SysMolecule state initialisation
		for(size_t i = 0; i < atom.size(); i++)
		{
			mol.atomxyz(i).setTo(0.0, 0.0, 0.0);
			atom[i].setRebuildRequired( true ); // mark atom as requiring attention by the Rebuilder
			atom[i].setKnownStructure( false ); // We dont know it until we have found an entry in the 'particles' array
		}

		// flag all file particles as unallocated
		std::vector<bool> allocatedParticle( particles.size(), false ); 

		bool changesOccured = true;
		// Loop over file particles
		while(changesOccured)
		{
			changesOccured = false;

			for( size_t i = 0; i < bioseq.size(); i++ )
			{
				int j = equiv[i];
				if( j != -1 )
				{
					const Sequence::ResidueInfo& ri = bioseq.getResidue(i);

					for(int k = mol.res[i].ifirst; k <= mol.res[i].ilast; k++)
					{
						ASSERT((atom[k].ir == i), CodeException,"Code failure in setPolymerPositions(), bioseq index doesnt match.");
						ASSERT((0 == ri.getFullName().compare(atom[k].parentname)),CodeException,"Code failure in setPolymerPositions(), bioseq resname doesnt match.");

						for(int m = particleStart[j]; m <= particleEnd[j]; m++)
						{
							if( allocatedParticle[m] || // check that we havent allocated this particle on a previous sweep
								0 != atom[k].pdbname.compare( particles[m].atomName ) // check that the atom names match!
								)
								continue;

							// If its not a backbone atom, some sidechain bonding needs to be checked!!
							// See problems described below!
							if( !atom[k].isBackbone() )
							{
								// Problem 1:
								// Now, in most cases, if the atom names are the same, we want to use the atom, BUT...
								// In the case of the example Leu -> Phe, both have CD1 and CD2 atoms, but Leus are tetrahedral
								// and Phes are planar. Therefore if the Leu coordinates are later used to rebuild the Phe, we will
								// end up with a distorted ring!! Therefore we only use this Fileparticle if the name matches 
								// AND the number of bonds is the same. This makes chemical sense...
								int molLookup = ffps.findMoleculeType( strucseq.getFullName(j) );
								ASSERT( molLookup != -1, ProcedureException, "Could not identify molecule type!!");
								const MoleculeDefinition& molParam = ffps.molecule[molLookup];
								int atomLookup = molParam.findAtomPDB( particles[m].atomName );
								if( atomLookup != -1 )
								{
									const AtomParameter& atomParam = molParam.atom[atomLookup];
									size_t lookupCovalency = atomParam.r_covalent.size();
									size_t covalency = atom[k].cov12atom.size();
									if( covalency != lookupCovalency && (covalency >=3 || lookupCovalency >= 3) )
									{
										if( verbose >= Verbosity::Normal )
										{
											printf("Ignorning due to covalency difference: " );
											particles[m].info();
										}
										// Its been assigned as processed, but it cant be as they dont match covalent types
										allocatedParticle[m] = true;
										changesOccured = true;
										continue;
									}
									else
									{
										// Problem 2:
										// Now sometimes as in the case of Tyr (or Phe) -> His one atom can be left "stranded"
										// and is therefore in a non-sensical location. In this example, I am refering to the CE1.
										// which in either structure is in a completely different relative position
										// to the remainder of the structure even though they are both ring structures!
										int builtPartnerCount = 0;
										for( size_t t1 = 0; t1 < atom[k].cov12atom.size(); t1++ )
										{
											int bondedPartnerIndex = atom[k].cov12atom[t1].i;
											bool partnerIsDefined = !atom[bondedPartnerIndex].isRebuildRequired();
											if(partnerIsDefined)
											{
												builtPartnerCount++;
											}
										}
										if( builtPartnerCount < 1 )
										{
											continue;
										}
									}
								}
							}

							if(!atom[k].isRebuildRequired())
							{
								// The atom has already been asigned by another 'particle'
								if(verbose)
								{
									printf("WARNING: Ambiguous particle equivalency for atom %i: ", k);
									atom[k].info();
								}
								continue; // use the 1st one that was encountered...
							}

							allocatedParticle[m] = true;
							changesOccured = true;

							//save particle's target position
							mol.atomxyz(k).setTo(particles[m].x, particles[m].y, particles[m].z);

							// set the relevent atom bit-flags
							atom[k].setRebuildRequired(false); // The Rebuilder can ignore this atom
							atom[k].setKnownStructure(true);
						}
					}
				}
			}
		}

		// Report what we have to build
		if(verbose)
		{
			for(size_t k = 0; k < atom.size(); k++)
			{
				if(atom[k].isRebuildRequired())
				{
					if( !atom[k].isHydrogen() )
					{
						printf("No structural data for heavy-atom: '%d'", k);
						atom[k].info();
					}
					else
					{
						// We would otherwise get a lot of info about missing 'H' atoms - we never really care -
						// missing heavy atoms is much more interesting output ... print only if we are screaming!
						if( verbose >= Verbosity::Scream )
						{
							printf("No structural data for hydrogen-atom: '%d'", k);
							atom[k].info();
						}
					}
				}
			}
		}

		return;
	}

	void setPolymerPositions( const FFParamSet& ffps, MoleculeBase& mol, const std::vector<FileParticle>& particles, Verbosity::Type verbose )
	{
		// setPolymerPositions() will cause each atoms pos states to be initialised either to the correct particle definition,
		// or to (0,0,0), with isKnownStructure() set to true and false respectively.

		// check we have something to work with...
		if( particles.size() == 0 || mol.nResidues() == 0 ) return;

		ParticleStore& atom = mol.atom;

		// SysMolecule state initialisation
		for(int i = 0; i < atom.size(); i++)
		{
			mol.atomxyz(i).setTo(0, 0, 0);
			atom[i].setRebuildRequired(true); // mark atom as requiring attention by the Rebuilder
			atom[i].setKnownStructure(false); // We dont know it until we have found an entry in the 'particles' array
		}

		const Sequence::BioSequence& seq = mol.getSequence();
		std::vector <bool> assigned(particles.size());  // remember the particles already assigned

		// scan through the list twice:
		//  round == 0:  compare the PDB names
		//  round == 1:  for unassigned particles/pdbatoms try matching the ffname (secondary name)
		for(int round = 0; round < 2; round++) 
		{
			int resSeqIndex = 0; // The 'particles' have no internal residue serial indexer, we have to count it here...
			char previousICode = particles[0].iCode;
			std::string previousName = particles[0].resName;
			int previousResNum = particles[0].resNum;

			// Attempt to assign each particle in turn to atoms in the SysMolcule
			for(int j = 0; j < particles.size(); j++)
			{
				// find atoms in the user-supplied PDB file
				if( (particles[j].iCode != previousICode) ||
					(particles[j].resNum != previousResNum)||
					(particles[j].resName.compare(previousName) != 0 ) )
				{
					previousICode = particles[j].iCode;
					previousResNum = particles[j].resNum;
					previousName = particles[j].resName;
					resSeqIndex++;
				}

				if( round == 0 )
				{	
					assigned[j] = false;
				}
				else
				{
					if( assigned[j] ) 
						continue; // on second round ignore unassinged atoms!
				}

				// Count through all the atoms looking for an equivelency. We also want to look for duplicate identities.
				for(int i = 0; i < atom.size(); i++)
				{
					if( round > 0 )
					{ 
						// on second round, ignore particles already having an assignement
						if(!atom[i].isRebuildRequired()) continue;
					}

					if( atom[i].ir != resSeqIndex )
						continue;
					if( 0 != seq.getResName(atom[i].ir).compare( particles[j].resName ) )
						continue;

					if( round == 0)
					{   
						// on round 1 compare to PDB names
						if( 0 != atom[i].pdbname.compare( particles[j].atomName ) )
							continue;
					}
					else
					{   
						// on round 2 compare to FF names
						if( 0 != atom[i].rawname.compare( particles[j].atomName ) )
							continue;
					}

					if(!atom[i].isRebuildRequired())
					{
						// The atom has already been asigned by another 'particle'
						if(verbose)
						{
							printf("WARNING: Ambiguous particle equivalency for atom %i\n", i);
						}
						continue; // use the 1st one that was encountered...
					}

					assigned[j] = true;

					//save particle's target position
					mol.atomxyz(i).setTo(particles[j].x, particles[j].y, particles[j].z);

					// set the relevent atom bit-flags
					atom[i].setRebuildRequired(false); // The Rebuilder can ignore this atom
					atom[i].setKnownStructure(true);
				}
			}
		}

		if(verbose)
		{
			for(int j = 0; j < particles.size(); j++)
			{
				if(!assigned[j])
				{
					printf("Unassigned particle '%d'", j);
					particles[j].info();
				}
			}
		}
	}


	// -----------------------------------------------
	// File_Molecule
	// -----------------------------------------------

	File_Molecule::File_Molecule( char _ChainID ) : chainID(_ChainID),
	Sequence::ExpSeqPair(bioSequence,structuralSequence) // assign the internal const references of ExpSeqPair to this classes member sequences
	{
	}

	File_Molecule::File_Molecule( const File_Molecule& _clone ) : Sequence::ExpSeqPair( bioSequence, structuralSequence)
	{
		bioSequence.append(_clone.bioSequence);
		structuralSequence.append(_clone.structuralSequence);
		atomLines = _clone.atomLines;
		chainID = _clone.chainID;

		// make sure we not only reassign the base classes internal sequence references to those belonging to this class
		// (above in the  Sequence::ExpSeqPair() constructor)
		// but also copy the current underlying alignment if there is one...
		CloneAlignment(_clone);
	}

	void File_Molecule::printToScreen() const
	{
		printf("ChainID: '%c'\n", chainID);
		if( bioSequence.size() == 0 )
		{
			printf("(Sequence Information is only available from the structural definition)\n");
			m_Seq2->printToScreen(' ',4,80);
			printf("\n\n");
		}
		else
		{
			// Call print on our base class, the alignment definition
			printf("(Sequence Information is available from SEQUENCE data)\n");
			Sequence::ExpSeqPair::printToScreen();
		}
	}

	void File_Molecule::setPrefix(Library::PrePostFix _First, Library::PrePostFix _Last)
	{
		structuralSequence.setPrefix(_First,_Last);
		bioSequence.setPrefix(_First,_Last);
	}


	FileBond::FileBond()
	{
		clear();
	}

	void FileBond::clear()
	{
		Type = Invalid;
		r1Name = "";
		r1ChainID = '\0';
		r1ResNum = -1;
		r1ICode = '\0';
		r2Name = "";
		r2ChainID = '\0';
		r2ResNum = -1;
		r2ICode = '\0';
	}

	bool FileBond::matches( const Sequence::BioSequence& _BioSeq, size_t _Index ) const
	{
		const Sequence::ResidueInfo& res = _BioSeq.getResidue(_Index);
		char iCode = res.getBiologicalICode();
		return
			( (
			0 == res.getFullName().compare(r1Name) &&
			_BioSeq.getChainID() == r1ChainID &&
			res.getBiologicalIndex() == r1ResNum &&
			(iCode == ' ' || iCode == r1ICode)
			) || (
			0 == res.getFullName().compare(r2Name) &&
			_BioSeq.getChainID() == r2ChainID &&
			res.getBiologicalIndex() == r2ResNum && // bioIndex is either the bioindex or the seqindex if bioIndex is not set
			(iCode == ' ' || iCode == r2ICode) // must be both as it is defaults to ' ' if it is not set
			) );
	}

	bool FileBond::matches( const FileParticle& atom ) const
	{
		return
			( (
			0 == atom.resName.compare(r1Name) &&
			atom.chainID == r1ChainID &&
			atom.resNum == r1ResNum &&
			atom.iCode == r1ICode
			) || (
			0 == atom.resName.compare(r2Name) &&
			atom.chainID == r2ChainID &&
			atom.resNum == r2ResNum &&
			atom.iCode == r2ICode
			) );
	}


	FileMoleculeMaker::FileMoleculeMaker( const FileImportBase& _parent )  : m_Parent(&_parent)
	{
	}

	FileMoleculeMaker::~FileMoleculeMaker()
	{
	}

	void FileMoleculeMaker::AssignBioSequence( File_Molecule& _Mol, const Sequence::BioSequence& _BioSequence ) const
	{
		_Mol.bioSequence.clear();
		if( _BioSequence.size() > 0 )
		{
			const FileCovalency& connections = m_Parent->getCovalency();
			if( connections.getCount() > 0 )
			{
				// check against the poxy-name dependent connections list to see if we have to rename our CYSs to CYXs (grrr)
				Sequence::ResidueInfo res(
					'\0',
					SIZE_T_FAIL,
					"");
				for( size_t j = 0; j < _BioSequence.size(); j++ )
				{
					res.setTo(_BioSequence.getResidue(j));
					for( size_t i = 0; i < connections.getCount(); i++ )
					{
						if( (connections[i].Type == FileBond::Disulphide) &&
							(connections[i].matches(_BioSequence,j)))
						{
							res.ForceRename("CYX");
							break;
						}
					}
					_Mol.bioSequence.append(res);
				}
			}
			else
			{
				// Just add the sequence as it is - isn't is so much simpler without 'FileCovalency'
				_Mol.bioSequence.append(_BioSequence);
			}
		}
	}

	bool FileMoleculeMaker::wouldAppendFileParticle( const FileParticle& _Particle ) const
	{
		static FileParticle part;
		part = _Particle;
		m_Parent->ReinterpretName( part ); // Virtual call to the parent: controls renaming of residues either from alias definitions or from some sort of knowledge
		return PassesFilter( part.resName );
	}

	bool FileMoleculeMaker::appendFileParticle( File_Molecule& _Mol, const FileParticle& _Particle ) const
	{
		static FileParticle part;
		part = _Particle;
		m_Parent->ReinterpretName( part ); // Virtual call to the parent: controls renaming of residues either from alias definitions or from some sort of knowledge
		if( PassesFilter( part.resName ) )
		{
			_Mol.atomLines.push_back(part);
			return true;
		}
		else
		{
			return false;
		}
	}

	void FileMoleculeMaker::DetectPrefix( File_Molecule& mol ) const
	{
		std::vector<std::string> foundClassNames; // which named classes do we have in our collection?
		std::vector<int> foundClassCounts; // count how many of each we can find.

		const Library::ClassMapper& classMap = getFileParent().getClass();
		const Sequence::BioSequence* seq = &mol.getSeq1();
		if( seq->size() == 0 )
		{
			// Then we are not using the Biological sequence - its empty
			// Use sequence2 - the structural sequence for type detection.
			seq = &mol.getSeq2();
		}

		if( seq->size() < 2 )
		{
			return; // no detection will be possible - we have a single molecule.
		}

		std::string className;
		for( size_t i = 0; i < seq->size(); i++ )
		{
			std::string resName = seq->getResName(i);		
			if( !classMap.getClassString(resName,className) ) continue;
			size_t foundAt = SIZE_T_FAIL;
			for( size_t j = 0; j < foundClassNames.size(); j++ )
			{
				if( 0 == className.compare(foundClassNames[j]))
				{
					foundAt = j;
					break;
				}
			}
			if( foundAt != SIZE_T_FAIL )
			{
				foundClassCounts[foundAt]++;
			}
			else
			{
				foundClassNames.push_back(className);
				foundClassCounts.push_back(1);
			}
		}

		if( foundClassCounts.size() == 0 )
		{
			printf("Warning: Failed to detect residue class, defaulting to NULL (no pre/postfix).\n");
			mol.setPrefix(Library::PrePostFix::NoFix,Library::PrePostFix::NoFix);
		}
		else if( foundClassCounts.size() == 1 )
		{
			printf("Polymer Type Detection Sucessful As: '%s'\n", foundClassNames[0].c_str() );
			Library::PrePostFix start;
			bool ok = classMap.getStartTerminusFromClass( foundClassNames[0], start );
			ASSERT(ok,CodeException,"ClassMapping oddness");
			Library::PrePostFix end;
			ok = classMap.getEndTerminusFromClass( foundClassNames[0], end );
			ASSERT(ok,CodeException,"ClassMapping oddness");
			mol.setPrefix(start,end);
		}
		else
		{
			// Ambiguous !!

			// Detect most common class.
			int max = -1;
			size_t maxAt = SIZE_T_FAIL;
			for( size_t j = 0; j < foundClassCounts.size(); j++ )
			{
				if( foundClassCounts[j] > max )
				{
					max = foundClassCounts[j];
					maxAt = j;
				}
			}
			printf("Warning: Failed to detect a sole residue class for molecule chain '%c'!\nPicking the most prevelant detection: '%s'\n", mol.getChainID(), foundClassNames[maxAt].c_str() );

			std::string startClass;
			if( !classMap.getClassString(seq->getFullName(0),startClass) ) 
			{
				printf("Warning: Start of polymer residue class is undetermined!\n");
			}
			else if( 0 != foundClassNames[maxAt].compare( startClass ) )
			{
				printf("Warning: Start of polymer residue class is not the same as the remainder of the polymer: '%s'!\n",startClass.c_str());
			}

			std::string endClass;
			if( classMap.getClassString(seq->getFullName(seq->size()-1),endClass) ) 
			{
				printf("Warning: End of polymer residue class is undetermined!\n");
			}
			else if( 0 == foundClassNames[maxAt].compare( startClass ) )
			{
				printf("Warning: End of polymer residue class is not the same as the remainder of the polymer: '%s'!\n",startClass.c_str());
			}

			Library::PrePostFix start;
			bool ok = classMap.getStartTerminusFromClass( foundClassNames[maxAt], start );
			ASSERT(ok,CodeException,"ClassMapping oddness");
			Library::PrePostFix end;
			ok = classMap.getEndTerminusFromClass( foundClassNames[maxAt], end );
			ASSERT(ok,CodeException,"ClassMapping oddness");
			mol.setPrefix(start,end);
		}
	}

	void FileMoleculeMaker::finaliseSequence( File_Molecule& mol ) const
	{
		std::vector<FileParticle>& atomLines = mol.atomLines;

		// Temporary Variables
		char previousICode = atomLines[0].iCode;
		std::string previousName = atomLines[0].resName;
		int previousResNum = atomLines[0].resNum;
		int sequentialIndex = 0;
		char singleID;

		// Here we use the alias mapper to obtain singleletter <-> long name mappings
		// NO renaming should or does occur in this function call ...
		const Library::AliasMapper& mapper = m_Parent->getAlias();

		mapper.lookupAlias(previousName,previousName);
		if( !mapper.lookupShortName(previousName,singleID) ) singleID = '?';
		Sequence::ResidueInfo res(singleID,sequentialIndex++,previousResNum,previousICode,previousName);

		for( size_t i = 0; i < atomLines.size(); i++ )
		{
			std::string resname = atomLines[i].resName;
			mapper.lookupAlias(resname,resname);

			if( (atomLines[i].iCode != previousICode) ||
				(atomLines[i].resNum != previousResNum)||
				(resname.compare(previousName) != 0 ) )
			{
				// add the res to the chain
				mol.structuralSequence.append( res );

				// residue has changed, change the 'previous' values
				previousICode = atomLines[i].iCode;
				previousResNum = atomLines[i].resNum;
				previousName = resname;

				// make the next one
				if( !mapper.lookupShortName(resname,singleID) ) singleID = '?';
				res = Sequence::ResidueInfo(singleID,sequentialIndex++,previousResNum,previousICode,previousName);
			}
		}

		// add the final residue
		mol.structuralSequence.append( res );

		// Detect the apropriate terminal residue prefixes for this polymer (e.g. N and C for a polypeptide)
		DetectPrefix(mol);

		// Perform alignment using disposable aligner
		m_Parent->getAligner().Align(mol);
		if( !mol.isValid() )
		{
			THROW(ProcedureException,"The sequence alignment produced from SEQRES and structural information is not consistent!\n");
		}
	}

	const FileImportBase& FileMoleculeMaker::getFileParent() const
	{
		return *m_Parent;
	}

	bool FileMoleculeMaker::isResequencableFilter() const
	{
		return m_Parent->isResequencableFilter();
	}

	bool FileMoleculeMaker::PassesFilter( const std::string& _ResName) const
	{
		return m_Parent->PassesFilter(_ResName);
	}

	std::string FileMoleculeMaker::getFilterText() const
	{
		return m_Parent->getFilterText();
	}





	FileImportBase::FileImportBase( const std::string& _filename, bool allocateClass )
	{
		m_Verbose = Verbosity::Silent;

		m_Filename = _filename;

		m_Covalency = &m_CovalencyDefault;
		m_Aligner = &m_AlignerDefault;

		m_FilterFlag = false;
		std::string m_FilterS = "";
		m_FilterE = Library::AllClasses;

		if( allocateClass ) m_Class = *Library::NamingConventions::getSingleton();
	}

	bool FileImportBase::ReinterpretName( FileParticle& atom ) const
	{
		ReinterpretFromCovalency( atom );
		return ReinterpretCore( atom );
	}

	void FileImportBase::ReinterpretFromCovalency( FileParticle& atom ) const
	{
		const FileCovalency& connections = getCovalency();
		for( size_t i = 0; i < connections.getCount(); i++ )
		{
			if( (connections[i].Type == FileBond::Disulphide) && (connections[i].matches(atom) ) )
			{
				// This needs to be changed to a CYX otherwise we will get
				// disulphides with CYS resides and an extra hydrogen atom knocking around
				atom.resName = "CYX";
				return;
			}
		}
	}

	bool FileImportBase::ReinterpretCore(FileParticle& atom ) const
	{
		// First allow the base class to use knowledge-based hacks if they are Enabled.
		// They are disabled by default...
		if( FileBabelBase::ReinterpretName(atom) ) return true;
		// Now use the current alias mapper to interpret the residue naming
		return getAlias().lookupAlias( atom.resName, atom.resName );
	}

	void FileImportBase::addAliasMapper(const Library::AliasMapper& _mapper)
	{
		m_Alias.addAliasMapper(_mapper);
	}

	void FileImportBase::setCovalency( FileCovalency& _covalency )
	{
		m_Covalency = &_covalency;
	}

	const std::string& FileImportBase::getFileName() const
	{
		return m_Filename;
	}

	Verbosity::Type FileImportBase::getVerbose() const
	{
		return m_Verbose;
	}

	void FileImportBase::setVerbosity( Verbosity::Type _BeVerbose )
	{
		m_Verbose = _BeVerbose;
	}

	const Sequence::AlignerBase& FileImportBase::getAligner() const
	{
		return *m_Aligner;
	}

	void FileImportBase::setAligner( Sequence::AlignerBase& _AlignerOverride )
	{
		m_Aligner = &_AlignerOverride;
	}

	void FileImportBase::setAlignerDefault()
	{
		m_Aligner = &m_AlignerDefault;
	}

	void FileImportBase::setAlignerDirect()
	{
		m_Aligner = &m_AlignerDirect;
	}

	void FileImportBase::setClassMapper(const Library::ClassMapper& _mapper)
	{
		// Copy over the current class definition
		m_Class = _mapper;
	}

	const Library::ClassMapper& FileImportBase::getClass() const
	{
		return m_Class;
	}

	Library::ClassMapper& FileImportBase::getClass()
	{
		return m_Class;
	}

	const Library::AliasMapperCollection& FileImportBase::getAlias() const
	{
		return m_Alias;
	}

	Library::AliasMapperCollection& FileImportBase::getAlias()
	{
		return m_Alias;
	}

	const FileCovalency& FileImportBase::getCovalency() const
	{
		return *m_Covalency;
	}

	FileCovalency& FileImportBase::getCovalency()
	{
		return *m_Covalency;
	}

	void FileImportBase::setFilter( const std::string& _filter )
	{
		m_FilterFlag = true;
		m_FilterS = _filter;
	}

	void FileImportBase::setFilter( Library::ResidueClass _filter )
	{
		m_FilterFlag = false;
		m_FilterE = _filter;
	}

	//void FileImportBase::setFilter( Library::ResidueClassGroup _filter )
	//{
	//	m_FilterFlag = false;
	//	m_FilterE = _filter;
	//}

	void FileImportBase::setFilter( Library::ExtdResidueClasses _filter )
	{
		m_FilterFlag = false;
		m_FilterE = _filter;
	}

	bool FileImportBase::isResequencableFilter() const
	{
		return !m_FilterFlag & isPolymerClass(m_FilterE);
	}

	bool FileImportBase::PassesFilter( const std::string& _ResName ) const
	{
		// We want to look at the class of the residue name that _ResName will resolve to,
		// NOT the cass that it currently is - peform a recursive lookup.
		std::string name = _ResName;
		std::string prevName;
		do
		{
			prevName = name;
			// Smile as the alias is recursively looked-up
			if( !getAlias().lookupAlias( name, name ) )
			{
				break;
			}
		}
		while( 0 != prevName.compare(name) );

		return m_FilterFlag ? getClass().isOfClass(m_FilterS,name) : getClass().isOfClass(m_FilterE,name);
	}

	std::string FileImportBase::getFilterText() const
	{
		return m_FilterFlag ? m_FilterS : getResidueClassString( m_FilterE );
	}




	FileInBase::FileInBase( const FFParamSet &_ffps, const std::string& _filename ) 
	: FileImportBase( _filename, false ), // false as we want to use the FFPS as our class source, see setClassMapper() in the constructor body!
	System( _ffps ),
	RebuilderOwner(),
	TreatMolsAsGroups(false),
	CentreOnBuild(false)
	{
		addAliasMapper(_ffps); // Always the **primary** alias mapper!
		addAliasMapper(*this); // This class IS a CustomAliasMapper, but importantly, still use the forcefield first.
		setClassMapper(_ffps);
	}

	bool FileInBase::ReinterpretName( FileParticle& atom ) const
	{
		ReinterpretFromCovalency( atom );

		if( -1 == ffps().findMoleculeType_withoutAlias(atom.resName) )
		{
			return FileImportBase::ReinterpretCore(atom);
		}
		else
		{
			// If the forcefield has a definition, then we never want to rename the residue, the forcefield definition
			// should always be used. e.g. We add a definition for pyrolysine to the forcefield, and a PDB file contains a
			// MODRES statement PLY --> LYS, then it is most desirable to use a proper pyro-lysine for simulations, rather
			// than "down-converting" to a LYS residue.
			return true;
		}
	}

	void FileInBase::SupplementaryBondScan( Molecule& sysMol )
	{
		// NOTE: This function is currently a bit HACK ey
		// This functionality is the 2nd phase of the SSBOND/FileCovalency support
		// Which supplies the CYS -> CYX conversions - however these should only
		// actually be bonded if in proximity. Maybee some further integration of this,
		// code is required, but that will be more difficult.... At the moment this works
		// well for proteins.

		// Proxies
		ParticleStore& atom = sysMol.atom;
		const Sequence::BioSequence& seq = sysMol.getSequence();

		// Find the indexes of CYX SG atoms that can form an SSBOND
		std::vector<size_t> cyxIndex;
		for( size_t j = 0; j < atom.size(); j++ )
		{		
			if( 0 == atom[j].pdbname.compare("SG") )
			{
				if( ( 0 == seq.getResName(atom[j].ir).compare("CYX") ) ||
					( 0 == seq.getResName(atom[j].ir).compare("NCYX") ) ||
					( 0 == seq.getResName(atom[j].ir).compare("CCYX") ) )
				{
					cyxIndex.push_back(j);
				}
			}
		}

		// Now bond the indexes in close proximity
		for( size_t i = 0; i < cyxIndex.size(); i++ )
		{
			for( size_t j = i + 1; j < cyxIndex.size(); j++ )
			{
				if( sysMol.atomxyz(cyxIndex[i]).sqrdist(sysMol.atomxyz(cyxIndex[j])) < 8.0 )
				{
					sysMol.addDisulphide(cyxIndex[i],cyxIndex[j]);
				}
			}
		}

		return;
	}

	void FileInBase::doCentring( Molecule& _NewMol )
	{
		Maths::dvector cog(0,0,0);
		int cnt = 0;
		for(int i = 0; i < _NewMol.nAtoms(); i++)
		{
			if( !_NewMol.atom[i].isRebuildRequired() )
			{
				cnt++;
				cog.add(_NewMol.atomxyz(i));
			}
		}
		if( cnt != 0 ) cog.div((double)cnt);
		for(int i = 0; i < _NewMol.nAtoms(); i++)
		{
			if( !_NewMol.atom[i].isRebuildRequired() )
			{
				_NewMol.atomxyz(i).sub(cog);
			}
			else
			{
				_NewMol.atomxyz(i).setTo(0,0,0);
			}
		}
	}

	void FileInBase::buildFinalise( Molecule& _NewMol, bool _IsPolymer )
	{
		// setPolymerPositions() will cause each atoms posRef() states to be initialised either to the correct particle definition,
		// or to (0,0,0), with isKnownStructure() set to true and false respectively.

		// Now, rebuild such things as missing sidechains and hydrogens using the current member builder - this is set to a default
		// value, but can be overidden from an external source.
		if( !getRebuilder().invokeBuild(_NewMol,getVerbose()) )
		{
			THROW(ProcedureException,"ERROR: Error occured during protein rebuild\n");
		}

		// Recentre the atoms we have found - Atoms requiring rebuild will not be considered.
		if( CentreOnBuild ) doCentring( _NewMol );

		// Missing atoms have a sensible defined position, either from structural data, or from a rebuild. This should now
		// be set as our true 'reference state'.
		_NewMol.setAsReference();

		// This is perhaps a bit of a hack, but I think this is the best way.
		// Basically all CYX residues should have been allocated in the 'File_Molecule' by now.
		// This *HAS* to be done before genpolymer(), but bonding *HAS* to occur afterwards.
		// In addition to this, CYX residues can be simulated without the bond (Mike said so...)
		// So they are only bonded here if they are in proximity to each other.
		if( _IsPolymer ) SupplementaryBondScan(_NewMol);

		/// add the resultant m_Molecule to the internal container of the 'System' base class
		add(_NewMol);
	}

	void FileInBase::build( File_Molecule& mol, bool _IsPolymer )
	{
		// 1) If we have a polymer, load the whole thing as a single molecule, containing multiple residues.
		// 2) If we have something like a solvent e.g. H20 we have 2 options: flagges by 'treatMolsAsGroups'
		//   mode a) Default: Each H20 "residue" is a single "molecule" in the eyes of the system.
		//   mode b) Non-standard and perhaps behaves unexpectedly: 
		//           Treat the solvent as a "non-bonded polymer" 
		//           i.e. All the H20 molecules belong to a and can be manipulated 
		//           as a single a "collection molecule"

		if(getVerbose()) printf("Building template molecule from structural data\n");

		if( _IsPolymer ) // Mode 1
		{		
			Molecule sysMol(ffps(),mol.getChainID());

			// make sure the coordinates reflect that of the available PDB file data, using the generated sequence allignment
			if(getVerbose()) printf("Assigning particles to built polymer using any available sequence alignment.\nMissing single atoms will be listed,entire missing sections \nwill not be part of an alignment.\n");
			if( mol.getBioLength() > 0 )
			{
				if( 0 != GeneratePolymer( mol.getBiologicalSeq(), sysMol, ffps(), getVerbose() ) != 0 ) 
				{
					throw ProcedureException("Error occured during polymer generation!");
				}
				// File_Molecule derives from 'Sequence::ExpSeqPair',
				// at this point the molecule 'ExpSeqPair' SHOULD have been aligned.
				ASSERT(mol.isSeqAligned(),CodeException,"Sequence must be aligned prior to FileInBase::build() call!");
				setPolymerPositions( ffps(), sysMol, mol, mol.getAtomLines(), getVerbose() ); // Use the complex version of 'setPolymerPositions' which uses a full sequence alignment for particle assignment
			}
			else
			{
				if( 0 != GeneratePolymer( mol.getStructuralSeq(), sysMol, ffps(), getVerbose() ) != 0 ) 
				{
					throw ProcedureException("Error occured during polymer generation!");
				}
				setPolymerPositions( ffps(), sysMol, mol.getAtomLines(), getVerbose() ); // Use the simple version of 'setPolymerPositions' which uses sequential indexing to assign particles
			}
			if(getVerbose()) printf("\n");

			buildFinalise( sysMol, true );
		}
		else if( TreatMolsAsGroups ) // Mode 2b
		{
			Molecule sysMol(ffps(),mol.getChainID());

			if( 0 != GenerateMoleculeSet( mol.getStructuralSeq(), sysMol, ffps(), getVerbose() ) )
			{
				throw ProcedureException("Error occured during molecule group generation!");
			}
			// make sure the coordinates reflect that of the available PDB file data
			if(getVerbose()) printf("setting to Coordinates from available particle data\n");
			setPolymerPositions( ffps(), sysMol, mol.getAtomLines(), getVerbose() ); // Use the simple version of 'setPolymerPositions' which uses sequential indexing to assign particles
			if(getVerbose()) printf("\n");

			buildFinalise( sysMol, false );
		}
		else // Mode 2a
		{
			// We must now re-interpret the file-atom list to chop it up into single molecules.

			const std::vector<FileParticle>& atomLines = mol.getAtomLines();
			size_t offset = 0; // offset from the start of the above atom list

			// Temporary Variables
			char previousICode = atomLines[0].iCode;
			std::string previousName = atomLines[0].resName;
			int previousResNum = atomLines[0].resNum;

			const Sequence::BioSequence& seq = mol.getStructuralSeq();
			for( size_t i = 0; i < seq.size(); i++ )
			{
				// For each "residue", make a brand new molecule
				Molecule newMol(ffps(),mol.getChainID());
				// Then find its particles in the full list
				std::vector<FileParticle> releventLines;

				if( 0 != loadMolecule(seq.getResidue(i),newMol,ffps()) )
				{
					throw ProcedureException("Error occured during molecule generation!");
				}

				// Find the required fileParticles
				for( /*none*/; offset < atomLines.size(); offset++ )
				{
					std::string resname = atomLines[i].resName;

					if( (atomLines[offset].iCode != previousICode) ||
						(atomLines[offset].resNum != previousResNum)||
						(resname.compare(previousName) != 0 ) )
					{
						// residue has changed, change the 'previous' values
						previousICode = atomLines[offset].iCode;
						previousResNum = atomLines[offset].resNum;
						previousName = resname;
						break; // break from the inner loop - we have found all lines for this molecule
					}
					else
					{
						releventLines.push_back(atomLines[offset]);
					}
				}

				setPolymerPositions( ffps(), newMol, releventLines, getVerbose() ); // Use the simple version of 'setPolymerPositions' which uses sequential indexing to assign particles

				buildFinalise( newMol, false );
			}
		}

		printf("Finished creating!\n\n");
	}

	void FileInBase::baseLoadCore( Library::ResidueClass _Class, FileMoleculeMaker& modelMaker, char chainID, const Sequence::BioSequence* _SequenceOverride )
	{	
		m_FilterE = _Class;
		File_Molecule mol = modelMaker.make(chainID,_SequenceOverride);

		if( !mol.isEmpty() )
		{
			if( Library::isPolymerClass( _Class ) )
			{
				mol.printToScreen(); // Only print the sequence if we have a polymer, otherwise we get things like "H20 H20 H20 .... H20"
				build( mol, true );
			}
			else
			{
				build( mol, false );
			}			
		}
	}

	void FileInBase::baseLoad( FileMoleculeMaker& modelMaker, char chainID, const Sequence::BioSequence* _SequenceOverride )
	{	
		// baseLoad() is an internal protected function call - we KNOW the chainID is valid.

		printf("Loading requested entities:\n");
		printf("---------------------------\n");

		if( m_FilterFlag ) // Check if string or enum filter-mode ?
		{
			// We are in filter-string mode - call single.
			File_Molecule mol = modelMaker.make(chainID,_SequenceOverride);
			mol.printToScreen();
			build( mol, true );
		}
		else
		{
			// We are in filter-enum mode
			// If flagged in the filter, load each present class in turn into a separate
			// molecule by enumerating the bit-flags
			Library::ExtdResidueClasses filter = m_FilterE; // Store the original filter
			int bitFlag = 1;
			for( int i = 0; i < Library::ResidueClassCount; i++ )
			{
				if( ((int)filter & bitFlag ) > 0 )   
				{
					baseLoadCore( (Library::ResidueClass) bitFlag, modelMaker, chainID, _SequenceOverride );
				}
				bitFlag *= 2; // Make it the next enum state... see Library::ResidueClasses defintion.
			}
			m_FilterE = filter; // Revert the old filter, it was changed prior to the above calls.
		}

		printf("\n");	
	}


	// -------------
	// FileToolBase
	// -------------

	FileToolBase::FileToolBase( const std::string& _filename ) : 
	FileImportBase(  _filename )
	{
	}


	// --------------
	// FileBabelBase
	// --------------

	FileBabelBase::FileBabelBase()
	{
		m_AtomNamer = Library::NamingConventions::getSingleton();

		m_UseKnowledge = false; // Knowledge-based hacking of names like 'MSE' and 'CSE' (and their atoms)

		// These 'Knowledge-based' terms are all disabled by default
		m_AtomNamer = Library::NamingConventions::getSingleton();
		m_AtomNameset = Library::PDB; // What convention would we be converting from?
		m_EnableAtomConversion = false; // Allow knowledge-based conversion of atom names from one convention to another?
	}

	bool FileBabelBase::ReinterpretName( FileParticle& atom ) const
	{
		if( m_UseKnowledge )
		{
			if( 0 == atom.resName.compare( "CSE" ) )
			{
				atom.resName = "CYS";
				if( 0 == atom.atomName.compare( "SE " ) )
					atom.atomName = " SD ";
				return true; // renaming has occured
			}
			else if( 0 == atom.resName.compare( "MSE" ) )
			{
				atom.resName = "MET";
				if( 0 == atom.atomName.compare( "SE " ) )
					atom.atomName = " SG ";
				return true; // renaming has occured
			}
		}
		if( m_EnableAtomConversion )
		{
			// coder warning - this function throws exceptions when it fails ...
			atom.atomName = m_AtomNamer->findAtomName(m_AtomNameset,atom.resName,atom.atomName);
			return true; // renaming has occured
		}
		return false;  // renaming hasn't occured
	}








	Load_PSF::Load_PSF( const FFParamSet &_ffps, const std::string& _filename ){


		std::vector < PSF_AtomLine > m_Atom;
		std::vector < IndexPair > m_Bonds;
		std::vector < IndexTriplet > m_Angles;

		enum Section { Pre, Atom, Bond, Angle, Torsion } current_section = Pre;

		// Scan file initially

		std::ifstream file(_filename.c_str(), std::ifstream::in);
		if( !file.is_open() ) throw(IOException( "Cannot open PSF file '" + _filename + "'!" ));

		std::string line;
		int line_count = 0;

		int psf_natoms;

		while(std::getline(file,line))
		{
			line_count ++;
			std::vector< std::string > token = chopstr( line, " \t" );
			printf(" %s ", line.c_str() );
			
			
			if( (token.size() >= 2 ) && ( !cmpstring( token[1] , "!NATOM" ) ) ){
				if( str2int( token[0], psf_natoms ) != 0){
					throw(IOException( "Number expected before !NATOM in PSF File: '" + _filename + "'!" ));		
				}
				current_section = Atom;
			}
			if( (token.size() >= 2 ) && ( !cmpstring( token[1] , "!NBOND" ) ) ){
				if( str2int( token[0], psf_natoms ) != 0){
					throw(IOException( "Number expected before !NBOND in PSF File: '" + _filename + "'!" ));		
				}
				current_section = Bond;
			}
			if( (token.size() >= 2 ) && ( !cmpstring( token[1] , "!NANGLE" ) ) ){
				if( str2int( token[0], psf_natoms ) != 0){
					throw(IOException( "Number expected before !NANGLE in PSF File: '" + _filename + "'!" ));		
				}
				current_section = Angle;
			}
			
			
			
			switch( current_section ){
				case Pre:
							if( line_count == 1 ){
								if( (token.size()==0) || ( !cmpstring( token[0] , "PSF" ) ) ){
									throw(IOException( "PSF file must begin with token 'PSF'. File: '" + _filename + "'!" ));
								}
							}



							break;
				case Atom:
							{
							PSF_AtomLine newatom;
							
							if( token.size() == 0 ) continue; // empty line

							if( token.size() < 9 ){
									throw(IOException( "Expected 9 Columns in PSF FILE: '" + _filename + "'!" ));		
									throw(ParseException( "Expected 9 Columns in Atom section '", _filename, line_count ));		
							}
							
							if(  str2int( token[0] , newatom.index ) != 0 ) throw(ParseException( "Expected integer in column 1 (atom index) ", _filename, line_count ) );
							newatom.segname = token[1];
							if(  str2int( token[2] , newatom.resnumber ) != 0 ) throw(ParseException( "Expected integer in column 3 (residue index) ", _filename, line_count ) );
							newatom.resname = token[3];
							newatom.atomname = token[4];
							newatom.ffname  = token[5];
							if(  str2double( token[6] , newatom.charge ) != 0 ) throw(ParseException( "Expected real number in column 7 (charge ) ", _filename, line_count ) );
							if(  str2double( token[7] , newatom.mass   ) != 0 ) throw(ParseException( "Expected real number in column 8 (mass ) ", _filename, line_count ) );
							if(  str2int  ( token[8] , newatom.extratoken) != 0 ) throw(ParseException( "Expected integer in column 9  ", _filename, line_count ) );
						
							m_Atom.push_back( newatom );
							}

							break;
						
				case Angle:
							break;
			}			


			





		}


		



		file.close(); // close the m_Stream to release the file_handle





		/// Now create the molecules !

		//Molecule newmolecule(sysspec.ffps()); // make ourselves a molecule to contain the imported system

		/*
		Sequence::BioSequence seq( sysspec.ffps() );
		int prevIR = INT_MAX;
		StringBuilder seqParse;
		*/
		
		/*
		for( size_t i = 0; i < m_Atom.size; i++ )
		{

		//		int iMol = sysspec.ffps().findMoleculeType( parentName );
		//		ASSERT(iMol!=-1,ProcedureException,"");
		//		int iAt = sysspec.ffps().molecule[iMol].findAtomPDB(pdbName);
		//		ASSERT(iAt!=-1,ProcedureException,"");
		//		const AtomParameter& atParam = sysspec.ffps().molecule[iMol].atom[iAt];
				
			Particle newparticle(atParam);

			newparticle.parentname = m_Atom.resname;

			if( parentName.size() > 3 ){
					newparticle.parentl3name = parentName.substr( parentName.size() - 3, 3 );
			}else{
					newparticle.parentl3name = parentName;
			}

			if( prevIR != sysDefs[i].parentnumber )
			{
				prevIR = sysDefs[i].parentnumber;
				seqParse.append( parentName );
				seqParse.append('-');
			}

			newparticle.parentl3name = parentName;
			newparticle.parentletter = parentName;

			newparticle.setKnownStructure(true);
			newparticle.setRebuildRequired(false);

				// residue number
				newparticle.ir = sysDefs[i].parentnumber;

				// read from our tra file
				traFile.read((char*)&getpos,sizeof(float)*3);
				newparticle.pos().setTo( getpos );

				newmolecule.addParticle(newparticle);
			}

			newmolecule.beginBonding();
			for( int i = 0; i < m_Header.atoms; i++ )
			{				
				for( int k = 0; k < sysDefs[i].n_cov12atoms; k++)
				{
					newmolecule.addBond( i, sysDefs[i].cov12atom[k] );
				}
			}
			newmolecule.endBonding();

			if( seqParse[seqParse.size()-1] == '-' )
			{
				seqParse.TruncateRightBy(1);
			}
			seq.setTo( seqParse.toString() );

			newmolecule.setBioSequence( seq );

			// and finally add the molecule to the parent container...
			sysspec.add(newmolecule);
		}

	*/



	}



} //namespace IO



