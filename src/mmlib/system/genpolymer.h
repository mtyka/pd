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

#ifndef __GENPOLYMER_H
#define __GENPOLYMER_H

namespace Sequence
{
	class BioSequence;
	class ResidueInfo;
}

class Molecule;
class MoleculeBase;
class FFParamSet;

#ifndef SWIG

/// For internal use only, but must be in the header as it is a friend of MoleculeBase
int GeneratorCore(
				  const Sequence::BioSequence &bioSeq,
				  Molecule &newmol,
				  const FFParamSet &ffps,
				  bool _Polymerise,
				  Verbosity::Type OutputLevel);

int GeneratePolymer(
					const std::string &aseq,              ///< Polymer Sequence
					Molecule &newmol,                     ///< Returns new molecule
					const FFParamSet &ffps,               ///< Forcefield Parameters
					Verbosity::Type OutputLevel);                         ///< make gurgling noises?

int GeneratePolymer(
					const Sequence::BioSequence &bioSeq,  ///< Polymer Sequence
					Molecule &newmol,                     ///< Returns new molecule
					const FFParamSet &ffps,               ///< Forcefield Parameters
					Verbosity::Type OutputLevel);                         ///< make gurgling noises?

int GenerateMoleculeSet(
						const std::string &aseq,              ///< Molecule-set Sequence
						Molecule &newmol,                     ///< Returns new molecule
						const FFParamSet &ffps,               ///< Forcefield Parameters
						Verbosity::Type OutputLevel);                         ///< make gurgling noises?

int GenerateMoleculeSet(
						const Sequence::BioSequence &bioSeq,  ///< Molecule-set Sequence
						Molecule &newmol,                     ///< Returns new molecule
						const FFParamSet &ffps,               ///< Forcefield Parameters
						Verbosity::Type OutputLevel);                         ///< make gurgling noises?

int loadMolecule(
				 const Sequence::ResidueInfo& _Res,
				 Molecule &newmol,
				 const FFParamSet &ffps);

/// Possibly DEPRECATED - this version does not assign the underlyings sequence == badness.
/// However there IS no sequence with a single molecule ...
int loadMolecule(
				 const std::string &MolName,
				 Molecule &newmol,
				 const FFParamSet &ffps);

#endif

Molecule PD_API NewProtein( FFParamSet &ffps, const std::string &aseq);

Molecule PD_API NewProteinHelix( FFParamSet &ffps, const std::string &aseq);

Molecule PD_API NewMolecule( FFParamSet &ffps, const std::string &MolName);

#endif


