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

#ifndef __RESIDUES_H
#define __RESIDUES_H

namespace Library
{
	// --------------------------------------------------------------------------------
	/// Enumeration of standard residue types and their grouped definitions
	/// Gratuitously 'borrowed' and enhanced from the UoB Framework ;-)
	// --------------------------------------------------------------------------------

	enum StandardResidues
	{
		Unknown = 0ul,

		// Protein
		A = 1ul, // Ala
		C = 2ul, // Cys
		D = 4ul, // Asp
		E = 8ul, // Glu
		F = 16ul, // Phe
		G = 32ul, // Gly
		H = 64ul, // His
		I = 128ul, // Ile
		K = 256ul, // Lys
		L = 512ul, // Leu
		M = 1024ul, // Met
		N = 2048ul, // Asn
		P = 4096ul, // Pro - Trans
		p = 8192ul, // Pro - Cis
		Q = 16384ul, // Gln
		R = 32768ul, // Arg
		S = 65536ul, // Ser
		T = 131072ul, // Thr
		V = 262144ul, // Val
		W = 524288ul, // Trp
		Y = 1048576ul, // Tyr

		// Nucleotide
		a = 2097152ul, // Adenine
		t = 4194304ul, // Thymine
		g = 8388608ul, // Guanine
		c = 16777216ul, // Cytosine
		u = 33554432ul // Uracil
	};

	// Coarse residue classes
	// Currently used in the sequence alignment classes
	enum StandardResidueProperties
	{
		AllProperties = A | C | D | E | F | G | H | I | K | L | M | N | P | p | Q | R | S | T | V | W | Y,
		NotGlyOrPro = A | C | D | E | F | H | I | K | L | M | N | Q | R | S | T | V | W | Y,

		// As defined on the ClustalW website:
		// http://www.ebi.ac.uk/clustalw/alignment_frame.html
		// Tiny = A | G | C | T | S,
		// Small = Tiny | N | D | V | P,
		// Aliphatic = V | I | L,
		// Aromatic = F | Y | W | H,
		// NonPolar = Aliphatic | Aromatic | M | C /*| K | H*/, // K and H ????
		// Positive = R | H | K,
		// Negative = D | E,
		// Charged = Positive | Negative,
		// Polar = Charged | T | S | C | N | Q | Y | W | H | K | R

		// Poster Definition
		// Small = G | A | S | T | C,
		// Hydrophobic = V | I | L | M | P,
		// Aromatic = F | Y | W,
		// Acidic = D | E,
		// Amide = N | Q,
		// Basic = H | K | R

		// My Older Definitions
		ShortHydrophobic = A | I | L | V,
		BulkyAromatic = F | H | W | Y,
		Polar = N | C | Q | S | T | Y,
		Positive = H | R | K, // H - debatable correspondence with K and R ... :-S ??
		Negative = D | E,
		Charged = Positive | Negative,

		Pyrimidine = t | c | u,
		Purine = a | g
	};

	// NOTE: These functions are case-sensitive!
	StandardResidues getResidue(char _QureyResidue);
	bool HasResidueProperty(StandardResidueProperties _Includes, StandardResidues _QureyResidue);
	bool HasResidueProperty(StandardResidueProperties _Includes, char _QureyResidue);
	bool HasSameResidueProperty( StandardResidues _QureyResidue1, StandardResidues _QureyResidue2 );
	bool HasSameResidueProperty( char _QureyResidue1, char _QureyResidue2 );
	bool HasSimilarResidueProperty( StandardResidues _QureyResidue1, StandardResidues _QureyResidue2 );
	bool HasSimilarResidueProperty( char _QureyResidue1, char _QureyResidue2 );
}

#endif

