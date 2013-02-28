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
#include "naming.h"
#include "tools/stringtool.h"

// ----------------------------------------------------------------------------------
// --- Convertion from aminoacid numbers to letter and vice versa ---
// ----------------------------------------------------------------------------------

const char *aalibfull[] = {
	"ALA","GLY","PRO","ARG",
	"ASN","ASP","CYS","GLN",
	"GLU","HIS","ILE","LEU",
	"LYS","MET","PHE","SER",
	"THR","TRP","TYR","VAL",

	"NALA","NGLY","NPRO","NARG",
	"NASN","NASP","NCYS","NGLN",
	"NGLU","NHIS","NILE","NLEU",
	"NLYS","NMET","NPHE","NSER",
	"NTHR","NTRP","NTYR","NVAL",

	"CALA","CGLY","CPRO","CARG",
	"CASN","CASP","CCYS","CGLN",
	"CGLU","CHIS","CILE","CLEU",
	"CLYS","CMET","CPHE","CSER",
	"CTHR","CTRP","CTYR","CVAL"
};

const char aalib[] = {
	'A','G','P','R',
	'N','D','C','Q',
	'E','H','I','L',
	'K','M','F','S',
	'T','W','Y','V',
	'A','G','P','R',
	'N','D','C','Q',
	'E','H','I','L',
	'K','M','F','S',
	'T','W','Y','V',
	'A','G','P','R',
	'N','D','C','Q',
	'E','H','I','L',
	'K','M','F','S',
	'T','W','Y','V'
};

//const char aalibnums[] = {
// 0,-1, 6, 5, 8,14,
// 1, 9,10,-1,12,11,
// 13, 4,-1, 2, 7,
// 3,15,16,-1,19,
// 17,-1,18
//};

char getAALetterFromFullName(char *resname){//This function
	//returns a one letter code from a 3 letter code
	for(int i = 0; i < 60; i++) {
		if(strcmp(aalibfull[i], resname) == 0)
			return aalib[i];
	}
	return '?';
}

char getAANumberFromFullName(char *resname){ //This function returns a number
	for(int i = 0; i < 60; i++) {
		if(strcmp(aalibfull[i], resname) == 0)
			return i;
	}
	return (char)-1;
}

char getAALetter(int aanum){
	if(aanum < 0)
		return '?';
	if(aanum >= 60)
		return '?';
	return (aalib[(aanum)]);
}

//const char *getAANameFull(int aanum){
// if(aanum < 0)
// return NULL;
// if(aanum >= 60)
// return NULL;
// return (aalibfull[(aanum)]);
//}

//int getAANumber(char aachar){
// if(aachar < 'A')
// return -1;
// if(aachar > 'Z')
// return -1;
// return (aalibnums[(aachar - 'A')]);
//}




// ----------------------------------------------------------------------------------
// --- Atom name convention converters ---
// ----------------------------------------------------------------------------------

const char Nameconvention_AAOrder[22] = { 'n', 'c',
'A', 'R', 'D', 'N',
'C', 'E', 'Q', 'G',
'H', 'I', 'L', 'K',
'M', 'F', 'P', 'S',
'T', 'W', 'Y', 'V'
};

char *PDB_Nameconvention[22][50] =
{
	// N terminus
	{"1H", "2H", "3H", "end"},
	// C terminus
	{"O", "OXT", "end"},
	// Ala, A, Alanine
	{"H", "HA", "1HB", "2HB", "3HB", "C", "CA", "CB", "N", "O", "end"},
	// Arg, R, Arginine
	{"H", "HA", "1HB", "2HB", "1HG", "2HG", "1HD", "2HD", "HE", "1HH1", "2HH1", "1HH2", "2HH2", "C", "CA", "CB", "CG", "CD", "CZ", "N", "NE", "NH1", "NH2", "O", "end"},
	// Asp, D, Aspartic Acid
	{"H", "HA", "1HB", "2HB", "C", "CA", "CB", "CG", "N", "O", "OD1", "OD2", "end"},
	// Asn, N, Aspargine
	{"H", "HA", "1HB", "2HB", "2HD2", "1HD2", "C", "CA", "CB", "CG", "N", "ND2", "O", "OD1", "end"},
	// Cys, C, Cysteine
	{"H", "HA", "1HB", "2HB", "HG", "C", "CA", "CB", "N", "O", "SG", "end"},
	// Glu, E, Glutamic Acid
	{"H", "HA", "1HB", "2HB", "1HG", "2HG", "C", "CA", "CB", "CG", "CD", "N", "O", "OE1", "OE2", "end"},
	// Gln, Q, Glutamine
	{"H", "HA", "1HB", "2HB", "1HG", "2HG", "2HE2", "1HE2", "C", "CA", "CB", "CG", "CD", "N", "NE2", "O", "OE1", "end"},
	// Gly, G, Glycine
	{"H", "1HA", "2HA", "C", "CA", "N", "O", "end"},
	// His, H, Histidine
	{"H", "HA", "1HB", "2HB", "HD1", "HD2", "HE1", "C", "CA", "CB", "CG", "CD2", "CE1", "N", "ND1", "NE2", "O", "end"},
	// Ile, I, isoleucine
	{"H", "HA", "HB", "1HG1", "2HG1", "1HG2", "2HG2", "3HG2", "1HD1", "2HD1", "3HD1", "C", "CA", "CB", "CG1", "CG2", "CD1", "N", "O", "end"},
	// Leu, L, Leucine
	{"H", "HA", "1HB", "2HB", "HG", "1HD1", "2HD1", "3HD1", "1HD2", "2HD2", "3HD2", "C", "CA", "CB", "CG", "CD1", "CD2", "N", "O", "end"},
	// Lys, K, Lysine
	{"H", "HA", "1HB", "2HB", "1HG", "2HG", "1HD", "2HD", "1HE", "2HE", "1HZ", "2HZ", "3HZ", "C", "CA", "CB", "CG", "CD", "CE", "N", "NZ", "O", "end"},
	// Met, M, Methionine
	{"H", "HA", "1HB", "2HB", "1HG", "2HG", "1HE", "2HE", "3HE", "C", "CA", "CB", "CG", "CE", "N", "O", "SD", "end"},
	// Phe, F, Phenyl Alanine
	{"H", "HA", "1HB", "2HB", "HD1", "HD2", "HE1", "HE2", "HZ", "C", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "N", "O", "end"},
	// Pro, P, Proline
	{"H2", "H1", "HA", "1HB", "2HB", "1HG", "2HG", "1HD", "2HD", "C", "CA", "CB", "CG", "CD", "N", "O", "end"},
	// Ser, S, Serine
	{"H", "HA", "1HB", "2HB", "HG", "C", "CA", "CB", "N", "O", "OG", "end"},
	// Thr, T, Threonine
	{"H", "HA", "HB", "HG1", "1HG2", "2HG2", "3HG2", "C", "CA", "CB", "CG2", "N", "O", "OG1", "end"},
	// Trp, W, Tryptophan
	{"H", "HA", "1HB", "2HB", "HD1", "HE1", "HE3", "HZ2", "HZ3", "HH2", "C", "CA", "CB", "CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2", "N", "NE1", "O", "end"},
	// Tyr, Y, Tyrosine
	{"H", "HA", "1HB", "2HB", "HD1", "HD2", "HE1", "HE2", "HH", "C", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "N", "O", "OH", "end"},
	// Val V, Valine
	{"H", "HA", "HB", "1HG1", "2HG1", "3HG1", "1HG2", "2HG2", "3HG2", "C", "CA", "CB", "CG1", "CG2", "N", "O", "end"}
};

char * XPLOR_Nameconvention[22][50] =
{
	// N terminus
	{"HT1", "HT2", "HT3", "end"},
	// C terminus
	{"OT1", "OT2", "end"},
	// Ala, A, Alanine
	{"HN", "HA", "HB1", "HB2", "HB3", "C", "CA", "CB", "N", "O", "end"},
	// Arg, R, Arginine
	{"HN", "HA", "HB2", "HB1", "HG2", "HG1", "HD2", "HD1", "HE", "HH11", "HH12", "HH21", "HH22", "C", "CA", "CB", "CG", "CD", "CZ", "N", "NE", "NH1", "NH2", "O", "end"},
	// Asp, D, Aspartic Acid
	{"HN", "HA", "HB2", "HB1", "C", "CA", "CB", "CG", "N", "O", "OD1", "OD2", "end"},
	// Asn, N, Aspargine
	{"HN", "HA", "HB2", "HB1", "HD21", "HD22", "C", "CA", "CB", "CG", "N", "ND2", "O", "OD1", "end"},
	// Cys, C, Cysteine
	{"HN", "HA", "HB2", "HB1", "HG", "C", "CA", "CB", "N", "O", "SG", "end"},
	// Glu, E, Glutamic Acid
	{"HN", "HA", "HB2", "HB1", "HG2", "HG1", "C", "CA", "CB", "CG", "CD", "N", "O", "OE1", "OE2", "end"},
	// Gln, Q, Glutamine
	{"HN", "HA", "HB2", "HB1", "HG2", "HG1", "HE21", "HE22", "C", "CA", "CB", "CG", "CD", "N", "NE2", "O", "OE1", "end"},
	// Gly, G, Glycine
	{"HN", "HA2", "HA1", "C", "CA", "N", "O", "end"},
	// His, H, Histidine
	{"HN", "HA", "HB2", "HB1", "HD1", "HD2", "HE1", "C", "CA", "CB", "CG", "CD2", "CE1", "N", "ND1", "NE2", "O", "end"},
	// Ile, I, isoleucine
	{"HN", "HA", "HB", "HG12", "HG11", "HG21", "HG22", "HG23", "HD11", "HD12", "HD13", "C", "CA", "CB", "CG1", "CG2", "CD1", "N", "O", "end"},
	// Leu, L, Leucine
	{"HN", "HA", "HB2", "HB1", "HG", "HD11", "HD12", "HD13", "HD21", "HD22", "HD23", "C", "CA", "CB", "CG", "CD1", "CD2", "N", "O", "end"},
	// Lys, K, Lysine
	{"HN", "HA", "HB2", "HB1", "HG2", "HG1", "HD2", "HD1", "HE2", "HE1", "HZ1", "HZ2", "HZ3", "C", "CA", "CB", "CG", "CD", "CE", "N", "NZ", "O", "end"},
	// Met, M, Methionine
	{"HN", "HA", "HB2", "HB1", "HG2", "HG1", "HE1", "HE2", "HE3", "C", "CA", "CB", "CG", "CE", "N", "O", "SD", "end"},
	// Phe, F, Phenyl Alanine
	{"HN", "HA", "HB2", "HB1", "HD1", "HD2", "HE1", "HE2", "HZ", "C", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "N", "O", "end"},
	// Pro, P, Proline
	{"HT2", "HT1", "HA", "HB2", "HB1", "HG2", "HG1", "HD2", "HD1", "C", "CA", "CB", "CG", "CD", "N", "O", "end"},
	// Ser, S, Serine
	{"HN", "HA", "HB2", "HB1", "HG", "C", "CA", "CB", "N", "O", "OG", "end"},
	// Thr, T, Threonine
	{"HN", "HA", "HB", "HG1", "HG21", "HG22", "HG23", "C", "CA", "CB", "CG2", "N", "O", "OG1", "end"},
	// Trp, W, Tryptophan
	{"HN", "HA", "HB2", "HB1", "HD1", "HE1", "HE3", "HZ2", "HZ3", "HH2", "C", "CA", "CB", "CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2", "N", "NE1", "O", "end"},
	// Tyr, Y, Tyrosine
	{"HN", "HA", "HB2", "HB1", "HD1", "HD2", "HE1", "HE2", "HH", "C", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "N", "O", "OH", "end"},
	// Val, V, Valine
	{"HN", "HA", "HB", "HG11", "HG12", "HG13", "HG21", "HG22", "HG23", "C", "CA", "CB", "CG1", "CG2", "N", "O", "end"}
};

char *PDAMBER_Nameconvention[22][50] =
{
	// N terminus
	{"HT1", "HT2", "HT3", "end"},
	// C terminus
	{"OT1", "OT2", "end"},
	// Ala, A, Alanine
	{"H", "HA", "HB1", "HB2", "HB3", "C", "CA", "CB", "N", "O", "end"},
	// Arg, R, Arginine
	{"H", "HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE", "HH11", "HH12", "HH21", "HH22", "C", "CA", "CB", "CG", "CD", "CZ", "N", "NE", "NH1", "NH2", "O", "end"},
	// Asp, D, Aspartic Acid
	{"H", "HA", "HB2", "HB3", "C", "CA", "CB", "CG", "N", "O", "OD1", "OD2", "end"},
	// Asn, N, Aspargine
	{"H", "HA", "HB2", "HB3", "HD21", "HD22", "C", "CA", "CB", "CG", "N", "ND2", "O", "OD1", "end"},
	// Cys, C, Cysteine
	{"H", "HA", "HB2", "HB3", "HG", "C", "CA", "CB", "N", "O", "SG", "end"},
	// Glu, D, Glutamic Acid
	{"H", "HA", "HB2", "HB3", "HG2", "HG3", "C", "CA", "CB", "CG", "CD", "N", "O", "OE1", "OE2", "end"},
	// Gln, Q, Glutamine
	{"H", "HA", "HB2", "HB3", "HG2", "HG3", "HE21", "HE22", "C", "CA", "CB", "CG", "CD", "N", "NE2", "O", "OE1", "end"},
	// Gly, G, Glycine
	{"H", "HA2", "HA3", "C", "CA", "N", "O", "end"},
	// His, H, Histidine
	{"H", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "C", "CA", "CB", "CG", "CD2", "CE1", "N", "ND1", "NE2", "O", "end"},
	// Ile, I, isoleucine
	{"H", "HA", "HB", "HG12", "HG13", "HG21", "HG22", "HG23", "HD11", "HD12", "HD13", "C", "CA", "CB", "CG1", "CG2", "CD1", "N", "O", "end"},
	// Leu, L, Leucine
	{"H", "HA", "HB2", "HB3", "HG", "HD11", "HD12", "HD13", "HD21", "HD22", "HD23", "C", "CA", "CB", "CG", "CD1", "CD2", "N", "O", "end"},
	// Lys, K, Lysine
	{"H", "HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "HE2", "HE3", "HZ1", "HZ2", "HZ3", "C", "CA", "CB", "CG", "CD", "CE", "N", "NZ", "O", "end"},
	// Met, M, Methionine
	{"H", "HA", "HB2", "HB3", "HG2", "HG3", "HE1", "HE2", "HE3", "C", "CA", "CB", "CG", "CE", "N", "O", "SD", "end"},
	// Phe, F, Phenyl Alanine
	{"H", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "HE2", "HZ", "C", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "N", "O", "end"},
	// Pro, P, Proline
	{"HN1", "HN2", "HA", "HB2", "HB3", "HG2", "HG3", "HD2", "HD3", "C", "CA", "CB", "CG", "CD", "N", "O", "end"},
	// Ser, S, Serine
	{"H", "HA", "HB2", "HB3", "HG", "C", "CA", "CB", "N", "O", "OG", "end"},
	// Thr, T, Threonine
	{"H", "HA", "HB", "HG1", "HG21", "HG22", "HG23", "C", "CA", "CB", "CG2", "N", "O", "OG1", "end"},
	// Trp, W, Tryptophan
	{"H", "HA", "HB2", "HB3", "HD1", "HE1", "HE3", "HZ2", "HZ3", "HH2", "C", "CA", "CB", "CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2", "N", "NE1", "O", "end"},
	// Tyr, Y, Tyrosine
	{"H", "HA", "HB2", "HB3", "HD1", "HD2", "HE1", "HE2", "HH", "C", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "N", "O", "OH", "end"},
	// Val, V, Valine
	{"H", "HA", "HB", "HG11", "HG12", "HG13", "HG21", "HG22", "HG23", "C", "CA", "CB", "CG1", "CG2", "N", "O", "end"}
};

int determineNamingConvention(){
	/*
	int *PDBconv;
	int *AMBERconv;
	int *XPLORconv;
	int PDBscore=0;
	int AMBERscore=0;
	int XPLORscore=0;
	int i;
	int *rtnum;
	int terminus;
	char rtletter;

	rtnum = new int[residues];
	PDBconv = new int[atoms];
	AMBERconv = new int[atoms];
	XPLORconv = new int[atoms];

	for(i=0;i<atoms;i++){
	terminus = 0;
	if(atom[i].resnum == firstres)

	rtletter = getAALetterFromFullName(&atom[i].resnam[0]);
	if(rtletter=='?') continue; // can't do anything
	for(i=2;i<22;i++){
	if(Nameconvention_AAOrder[i] == rtletter){ rtnum = i; break; }
	}
	if(i==22){
	printf("ERROR: Can't find PDB Naming convention conversion table \n");
	continue;
	}

	// PDB
	PDBconv[i] = AtomName_FindInPDB(rtnum,&atom[i].atnam[0]);
	if(PDBconv[i]>=0) PDBscore++;

	// AMBER
	AMBERconv[i] = AtomName_FindInAMBER(rtnum,&atom[i].atnam[0]);
	if(AMBERconv[i]>=0) AMBERscore++;

	// AMBER
	XPLORconv[i] = AtomName_FindInXPLOR(rtnum,&atom[i].atnam[0]);
	if(XPLORconv[i]>=0) XPLORscore++;
	}

	if( (PDBscore > AMBERscore) && (PDBscore > XPLORscore) ){


	}else
	if( (AMBERscore > PDBscore) && (AMBERscore > XPLORscore) ){


	}else
	if( (XPLORscore > PDBscore) && (XPLORscore > AMBERscore) ){


	}else



	return -1; // if we get here, we didn't find it
	}

	*/
	return 0;
}

int AtomName_FindInPDB(int accessnr, char *iname){
	char * newname;
	int j;
	for(j = 0; j < 50; j++) {
		newname = PDB_Nameconvention[accessnr][j];
		if(strcmp(newname, "end") == 0) {
			newname = NULL;
			return -1;
		}
		if(strcpy(newname, iname) == 0)
			return j;
	}
	return -1;
}

int AtomName_FindInXPLOR(int accessnr, char *iname){
	char * newname;
	int j;
	for(j = 0; j < 50; j++) {
		newname = PDB_Nameconvention[accessnr][j];
		if(strcmp(newname, "end") == 0) {
			newname = NULL;
			return -1;
		}
		if(strcpy(newname, iname) == 0)
			return j;
	}
	return -1;
}

int AtomName_FindInAMBER(int accessnr, char *iname){
	char * newname;
	int j;
	for(j = 0; j < 50; j++) {
		newname = PDB_Nameconvention[accessnr][j];
		if(strcmp(newname, "end") == 0) {
			newname = NULL;
			return -1;
		}
		if(strcpy(newname, iname) == 0)
			return j;
	}
	return -1;
}

int AtomName_PDAMBER2PDB(char *resname, char *iname, char *oname){
	char *oldname = NULL;
	char *newname = NULL;
	int i, j;
	int accessnr;

	if(strcmp(resname, "CTER") == 0) {
		accessnr = 0;
	} else if(strcmp(resname, "NTER") == 0) {
		accessnr = 1;
	} else {
		char rtletter;
		rtletter = getAALetterFromFullName(resname);
		if(rtletter == '?') {
			strcpy(oname, iname); // just keep name - can't do anything
			return -1;
		}
		//printf("rtletter is '%c' \n",rtletter);
		for(i = 2; i < 22; i++) {
			if(Nameconvention_AAOrder[i] == rtletter) {
				accessnr = i;
				break;
			}
		}
		if(i == 22) {
			printf("ERROR: Can't find PDB Naming convention conversion table \n");
			strcpy(oname, iname); // just keep name - can't do anything
			return -1;
		}
	}

	for(i = 0; i < 50; i++) {
		oldname = PDAMBER_Nameconvention[accessnr][i];
		if(strcmp(oldname, "end") == 0) {
			oldname = NULL;
			break;
		}
		if(strcmp(oldname, iname) == 0)
			break;
	}

	if(oldname == NULL) { // didn't find atom name
		strcpy(oname, iname); // just keep name - can't do anything
		return -1;
	}

	for(j = 0; j < 50; j++) {
		newname = PDB_Nameconvention[accessnr][j];
		if(strcmp(newname, "end") == 0) {
			newname = NULL;
			break;
		}
		if(j == i) { // found it !
			strcpy(oname, newname);
			return 0;
		}
	}

	strcpy(oname, iname); // just keep name - can't do anything
	return -1; // if we get here, we didn't find it
}

