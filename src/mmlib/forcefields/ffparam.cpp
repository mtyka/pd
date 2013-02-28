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

// Pre-compiled
#include "global.h"

// Self
#include "ffparam.h"

// MMLib
#include "tools/stringbuilder.h"
#include "tools/io.h"
#include "library/libpath.h"
#include "library/elements.h"
#include "forcefield.h"

using namespace Maths;
using namespace Physics;
using namespace Library;

// The global deliminator for the forcefield format
const char *wspace = " \t\12\15";
const char *commentchar = "#\12\15";

const FFParamSet nullffps;
const MoleculeDefinition nullmolecule;
const AtomParameter nullatomparameter;

// takes identifiers like "-CA" and splits them into -1 (residue offset) and "CA" (atomname)
// currently a simple version allowing only - or + as prefixes, maybe we'll extend this
// one day to a more general treatment
int interpretAtomReference(const std::string &identifier,
                           std::string &atomname,
                           int &residueoffset)
{
    if(identifier[0] == '-')
    {
        atomname = identifier.substr(1); // miss out 1st character
        residueoffset = (int) -1;
    }
    else if(identifier[0] == '+')
    {
        atomname = identifier.substr(1); // miss out 1st character
        residueoffset = (int) +1;
    }
    else
    {
        atomname = identifier; // miss out 1st character
        residueoffset = (int) 0;
    }

    return 0;
}

/// Asserts that a given token contains ONLY alpha numeric characters
/// NOTE: That names, particularily residue names may not contain either
/// * ( ) or - (due to string formatting opperations in biocseqeunce.)
/// However _ IS allowed in residue names !!
bool token2name(const std::string& token, std::string& name )
{
    for( size_t i = 0; i< token.size(); i++ ) 
        if( !isalnum(token[i]) && (token[i]!='_') ) 
            return false;	
    name = token;
    return true;
}

/// Asserts that a given token contains ONLY allowed characters
bool token2name(const std::string& token, std::string& name, const std::string& _forbiddenChars )
{
    for( size_t i = 0; i< token.size(); i++ ) 
        for( size_t j = 0; j< token.size(); j++ ) 
            if( token[i] == _forbiddenChars[j] ) 
                return false;	
    name = token;
    return true;
}

/// --------------------------------------------------------------
// BondTypeParameter Class Function Definitions
/// reads out a line like this and saves it in itself
/// BOND BB_C BB_O 1.2330 9.2000
int BondTypeParameter::readDefinitionLine(const std::string linestring, FFParamSet & ffps, ParseErrorLogger &errlog){

    std::vector<std::string> token;
    token = chopstr(linestring,wspace);

    if(token.size() < 5)
    {
        errlog.logError("Insufficient BOND parameters - expected two atom type paramters, followed by bondlength, forceconstant and bondorder");
        return -1;
    }

    used = 1; // set used flag to 'yes'
    i = ffps.findAtomType(token[1]);
    if(i < -1)
    {
        errlog.logError("Type identifier '" + token[1] + "' unknown ");
        return -1;
    }
    j = ffps.findAtomType(token[2]);
    if(j < -1)
    {
        errlog.logError("Type identifier '" + token[2] + "' unknown ");
        return -1;
    }

    if(str2double(token[3],length)!=0)
    { errlog.logError("SYNTAX ERROR: Float expected for bond length"); return -1; }
    if(str2double(token[4],forceconstant)!=0)
    { errlog.logError("SYNTAX ERROR: Float expected for bond forceconstant"); return -1; }
    if(str2double(token[5],bondorder)!=0)
    { errlog.logError("SYNTAX ERROR: Float expected for bond bondorder"); return -1; }

    //Do any unit conversions
    forceconstant *= Physics::PhysicsConst::kcal2J / PhysicsConst::Na; // kcal/mol to J (SI)

    return 0;
}

/// compares if the bond Type matches t{i,j}
int BondTypeParameter::cmpBond(int ti, int tj) const
{ 
    int result;
    int wildcards;

    result = 2;
    wildcards = 0;

    if(i == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(i == ti)
        result -= 1; // otherwise look for a match
    if(j == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(j == tj)
        result -= 1; // otherwise look for a match

    if(result == 0) // if all 2 matched, then result is 0 by now
        return wildcards; // return number of wildcards

    // if no forward match is achieved, try backwards

    result = 2;
    wildcards = 0;

    if(i == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(i == tj)
        result -= 1; // otherwise look for a match
    if(j == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(j == ti)
        result -= 1; // otherwise look for a match

    if(result == 0) // if all 2 matched, then result is 0 by now
        return wildcards; // return number of wildcards

    return -1; // otherwise return -1 indicating no match
}



/// --------------------------------------------------------------
// AngleTypeParameter Class Function Definitions

int AngleTypeParameter::readDefinitionLine(const std::string linestring, FFParamSet & ffps, ParseErrorLogger &errlog){
    // reads out a line like this and saves it in itself
    // ANGLE BB_N BB_C BB_O 1.2330 9.2000

    std::vector<std::string> token;
    token = chopstr(linestring,wspace);

    if(token.size() < 6)
    {
        errlog.logError("Insufficient ANGLE parameters - expected three atom type identifiers, followed by equilibrium angle and forceconstant");
        return -1;
    }

    used = 1; // set used flag to 'yes'
    i = ffps.findAtomType(token[1]);
    if(i < -1)
    {
        errlog.logError("Type identifier  '" + token[1] + "' unknown ");
        return -1;
    }
    a = ffps.findAtomType(token[2]);
    if(a < -1)
    {
        errlog.logError("Type identifier  '" + token[2] + "' unknown ");
        return -1;
    }
    j = ffps.findAtomType(token[3]);
    if(j < -1)
    {
        errlog.logError("Type identifier  '" + token[3] + "' unknown ");
        return -1;
    }

    if(str2double(token[4],angle)!=0)
    { errlog.logError("SYNTAX ERROR: Float expected for angle"); return -1; }
    if(str2double(token[5],forceconstant)!=0)
    { errlog.logError("SYNTAX ERROR: Float expected for angle forceconstant"); return -1; }

    //Do any unit conversions
    angle *= Maths::MathConst::PIOver180; // degrees to radians
    forceconstant *= PhysicsConst::kcal2J / PhysicsConst::Na; // kcal/mol/rad to J/rad (SI)

    return 0;
}

int AngleTypeParameter::cmpAngle(int ti, int ta, int tj) const{ // compares if the angle Type matches t{i,a,j}
    int result;
    int wildcards;

    result = 3;
    wildcards = 0;

    if(i == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(i == ti)
        result -= 1; // otherwise look for a match
    if(a == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(a == ta)
        result -= 1; // otherwise look for a match
    if(j == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(j == tj)
        result -= 1; // otherwise look for a match

    if(result == 0) // if all 4 matched, then result is 0 by now
        return wildcards; // return number of wildcards

    // if no forward match is achieved, try backwards

    result = 3;
    wildcards = 0;

    if(i == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(i == tj)
        result -= 1; // otherwise look for a match
    if(a == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(a == ta)
        result -= 1; // otherwise look for a match
    if(j == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(j == ti)
        result -= 1; // otherwise look for a match

    if(result == 0) // if all 4 matched, then result is 0 by now
        return wildcards; // return number of wildcards

    return -1; // otherwise return -1 indicating no match
}


/// --------------------------------------------------------------
// TorsionTypeParameter Class Function Definitions

// compares if the torsion Type matches t{i,a,b,j}
int TorsionTypeParameter::cmpTorsion(int ti, int ta, int tb, int tj, int ttype) const
{
    //first check if Type is equal (don't compare torsions with impropers)
    if(Type != ttype)
        return -1;

    //otherwise check correspondance
    int result;
    int wildcards;

    result = 4;
    wildcards = 0;

    if(i == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(i == ti)
        result -= 1; // otherwise look for a match
    if(a == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(a == ta)
        result -= 1; // otherwise look for a match
    if(b == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(b == tb)
        result -= 1; // otherwise look for a match
    if(j == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(j == tj)
        result -= 1; // otherwise look for a match

    if(result == 0) // if all 4 matched, then result is 0 by now
        return wildcards; // return number of wildcards

    // if no forward match is achieved, try backwards

    result = 4;
    wildcards = 0;

    if(i == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(i == tj)
        result -= 1; // otherwise look for a match
    if(a == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(a == tb)
        result -= 1; // otherwise look for a match
    if(b == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(b == ta)
        result -= 1; // otherwise look for a match
    if(j == -1)
    {
        result -= 1;
        wildcards++;
    } // if its a wildcard, its always a match
    else if(j == ti)
        result -= 1; // otherwise look for a match

    if(result == 0) // if all 4 matched, then result is 0 by now
        return wildcards; // return number of wildcards

    return -1; // otherwise return -1 indicating no match
}

void TorsionTypeParameter::zero()
{
    terms = 0;
    a = -2;
    i = -2;
    j = -2;
    b = -2;
}

int TorsionTypeParameter::addTerm(double newVn, double newn, double newgamma)
{
    // Look if a term with this multiplicity has already been added
    int newterm = terms;

    for( int iterm = 0; iterm < terms; iterm ++ )
    {
        if( n[iterm] == newn ){
            newterm = iterm;
            info();
            printf("INFO: Patching this torsion with new term with multiplicity n = %.0f\n", newn);
            break;	
        }
    }

    if(terms >= 4)
    {
        info();
        throw(ParseException("Maximum of 4 fourier terms supported per given torsion"));
    }

    Vn[newterm] = newVn;
    n[newterm] = newn;
    gamma[newterm] = newgamma;

    //printf("G: %f  %f  %f  \n", Vn[newterm], n[newterm],  gamma[newterm] ); 
    if( newterm == terms ) terms++;

    return 0;
}

/// reads out a line like this and adds it to itself
/// TORSION BB_CA BB_N BB_C BB_O ??? ??? ???
///
/// list is an array, of length list_length
int TorsionTypeParameter::readDefinitionLine(
    std::vector<TorsionTypeParameter> &list,
    int list_length,
    const std::string &linestring,
    FFParamSet &ffps,
    ParseErrorLogger &errlog
    )
{

    int prevtors;
    double newVn, newn, newgamma;

    std::vector<std::string> token;
    token = chopstr(linestring,wspace);

    if(token.size() < 8) {
        errlog.logError("Insufficient TORSION parameters - expected four atom identifiers, followed by the scaling parameterV, the offset gamma parameter and the multiplicity n");
        return -1;
    }

    //determine if it's a torsion or an improper dihedral
    if(cmpstring(token[0], "TORSION"))
    {
        Type = 0;
    } else if(cmpstring(token[0], "IMPROPER"))
    {
        Type = 1;
    }
    else
    {
        THROW(CodeException,"CODE ERROR: TorsionTypeParameter::readDefinitionLine(...):: neither IMPROPER nor TORSION ?? ");
    }

    used = 1; // set used flag to 'yes'
    i = ffps.findAtomType(token[1]);
    if(i < -1)
    {
        errlog.logError("Type identifier '" + token[1]  + "' unknown"); 
        return -1;
    }              
    a = ffps.findAtomType(token[2]);                     
    if(a < -1)
    {   
        errlog.logError("Type identifier '" + token[2]  + "' unknown");
        return -1;
    }              
    b = ffps.findAtomType(token[3]);                     
    if(b < -1)
    {   
        errlog.logError("Type identifier '" + token[3]  + "' unknown");
        return -1;
    }              
    j = ffps.findAtomType(token[4]);                     
    if(j < -1)
    {   
        errlog.logError("Type identifier '" + token[4]  + "' unknown"); 
        return -1;
    }

    if(str2double(token[5],newVn)!=0)
    { 
        errlog.logError("SYNTAX ERROR: Float expected for Vn");
        return -1; 
    }
    if(str2double(token[6],newgamma)!=0)
    { 
        errlog.logError("SYNTAX ERROR: Float expected for Gamma");
        return -1; 
    }
    if(str2double(token[7],newn)!=0)
    { 
        errlog.logError("SYNTAX ERROR: Float expected for n");
        return -1;
    }

    // Do any unit conversion
    newVn *= PhysicsConst::kcal2J / PhysicsConst::Na; // kcal/mol --> J / atom
    newgamma *= Maths::MathConst::PIOver180; //deg --> rad

    //now find out if there was a previous torsion of the same Type (i.e.
    //we need to add a new fourier term, not an entirely new torsion Type
    for(prevtors = 0; prevtors < list_length; prevtors++)
    {
        if( (list[prevtors].i == i) &&
            (list[prevtors].a == a) && 
            (list[prevtors].b == b) && 
            (list[prevtors].j == j) && 
            (list[prevtors].Type == Type)) 
        {
            break;
        }
    }

    if(prevtors < list_length)
    {
        // yep there was a previous torsion of same Type
        list[prevtors].addTerm(newVn, newn, newgamma); // ok, add to previous one
        return 1; // signal not increase torsion counter
    }
    else
    {
        // no previous torsion of same Type found
        addTerm(newVn, newn, newgamma);
        return 0; // signal to increase torsion counter
    }
}

void TorsionTypeParameter::info()
{
    int h;

    printf("TORSION %3d %3d %3d %3d %6.3lf %6.1lf %2.lf\n",
        i, a, b, j, Vn[0] * PhysicsConst::J2kcal * PhysicsConst::Na, gamma[0] * Maths::MathConst::OneEightyOverPI, n[0]);
    for(h = 1; h < terms; h++)
    {
        printf(" %6.3lf %6.1lf %2.lf\n", Vn[h] * PhysicsConst::J2kcal * PhysicsConst::Na, gamma[h] * Maths::MathConst::OneEightyOverPI, n[h]);
    }
}


/// --------------------------------------------------------------
// MoleculeDefinition Class Function Definitions

int MoleculeDefinition::findAtomRaw(const std::string &atom_name) const
{
    int index = -1; // sift through all atoms of that aacid to find the one referenced
    for(size_t i = 0; i < atom.size(); i++)
    {
        if(cmpstring(atom[i].rawname, atom_name))
        {
            index = i;
            break;
        }
    }
    return index;
}

int MoleculeDefinition::findAtomPDB(const std::string &atom_name) const
{
    int index = -1; // sift through all atoms of that aacid to find the one referenced
    for(size_t i = 0; i < atom.size(); i++)
    {
        if(cmpstring(atom[i].pdbname, atom_name))
        {
            index = i;
            break;
        }
    }
    return index;
}


int MoleculeDefinition::decodeRestraintList(
    const std::string &data,
    std::vector<CovalentLink> &bondlist,
    ParseErrorLogger &errlog 
    )
{
    std::vector<std::string> restlist;
    restlist = chopstr(data,","); // Chop up string by deliminator

    bondlist.clear();

    unsigned i;
    for(i=0;i<restlist.size();i++)
    {
        CovalentLink newbond;

        // interpret atom identifier (disect out preceeding - or +)
        if(interpretAtomReference(restlist[i], newbond.ani, newbond.roi))
        {
            errlog.logError("SYNTAX ERROR: Dont understand " + restlist[i] );
            return -1;
        }

        // Now find internal index (adapted OLD code), -1 if the identifier
        // refers to the next/previous residue

        // First check if the later is true by checking the residue offset of i
        // if it's not 0 the atom name refers to an atom in a different residue,
        // hence internal atom numbers don't make sense: set it to -1 to indicate that
        // and move on to next bond

        if(newbond.roi != 0)
        {
            newbond.i = -1;

            // Jon Removed - frwdlink/backlink should be defined by fffile ?!?
            // Set link properties of molecule
            //if(newbond.roi < 0) backlink = true; // at least 1 back link is present
            //if(newbond.roi > 0) frwdlink = true; // at least 1 forward link is present
        }
        else
        {
            int rest_index = findAtomRaw(newbond.ani);
            if(rest_index == -1)
            { 
                // If atom name unknown, warn & skip
                char temp[2] = { letter, 0 };
                errlog.logError("ERROR: Unknown atom identifier in molecule/residues " + std::string(&temp[0]) + "/" + l3_name + "/" + name + " : " + newbond.ani );
                continue;
            }
            newbond.i = rest_index;
        }

        bondlist.push_back(newbond);
    }

    return 0;
}


int MoleculeDefinition::checkRestraintList(  ParseErrorLogger &errlog )
{
    unsigned j, k, l;
    int resti, ok;
    int return_value = 0;

    for(j = 0; j < atom.size(); j++)
    {
        // for each atom in that amino acid
        for(k = 0; k < atom[j].r_covalent.size(); k++)
        {
            // for every covalent restraint

            // if it's not 0 the atom name refers to an atom in a different residue,
            // hence we cant check the back link at this stage, this will have to be done
            // later once the entire molecule is assembled
            if(atom[j].r_covalent[k].roi != 0) continue;

            // get the atom reference number for this restraint
            resti = atom[j].r_covalent[k].i;

            // Now check if the atom linked to, also links back to this atom.
            ok = 0; // (false)
            for(l = 0; l < atom[resti].r_covalent.size(); l++)
            {
                if(atom[resti].r_covalent[l].i == j)
                {
                    // Is atom j part of atom rest i's restraint list? Yes?
                    ok = 1; // (true)
                    break; // break outermost for loop (l)
                }
            }

            if(ok == 0)
            {
                errlog.logError("ERROR: Missing back link for restraint " + int2str(k + 1) + ", atom " + atom[j].pdbname + ", molecule " + std::string( &name[0] ) );
                return_value = -1; // error found
            }
        }
    }
    return return_value;
}

int MoleculeDefinition::decodeRestraintList( ParseErrorLogger &errlog )
{
    double totalcharge = 0;
    for(size_t atnum = 0; atnum < atom.size(); atnum++)
    {
        decodeRestraintList(atom[atnum].restraint, atom[atnum].r_covalent, errlog);
        totalcharge += atom[atnum].charge;
    }

    // check that total charge is integral
    double chargedev1 = fabs(floor(totalcharge) - totalcharge);
    double chargedev2 = fabs(ceil(totalcharge) - totalcharge);
    double chargedev = min(chargedev1, chargedev2);

    if(chargedev > 0.05)
    {
        errlog.logError("WARNING: Molecule/Residue " + name + " has non-integral total charge " + double2str( totalcharge ) );
    }

    // Now do a check of the restraint information looking for errors
    // in the definition file
    if(checkRestraintList( errlog ) != 0)
    {
        errlog.logError("ERROR: Errors found in restraint lists - check definition files \n");
        return -1;
    }

    return 0;
}


/// reads out a line like this and adds it to itself
/// TORSION/IMPROPER +CA -N C O
int DihedralDefinition::readDefinitionLine(const std::string &linestring, ParseErrorLogger &errlog)
{
    std::vector<std::string> token;
    token = chopstr(linestring,wspace);

    if(token.size() < 5)
    {
        errlog.logError("Insufficient DIHEDRAL parameters (" + int2str(token.size()) + " given) "); 
        return -1;
    }

    //determine if it's a torsion or an improper dihedral
    if(cmpstring(token[0], "TORSION") || cmpstring(token[0], "DELTORSION"))
    {
        Type = 0;
    } 
    else if(cmpstring(token[0], "IMPROPER") || cmpstring(token[0], "DELIMPROPER"))
    {
        Type = 1;
    } 
    else 
    {
        THROW(CodeException,"CODE ERROR: TorsionTypeParameter::readDefinitionLine(...):: neither IMPROPER nor TORSION ?? ");
        return -1;
    }

    if(interpretAtomReference(token[1], ani, roi) != 0) 
    {
        errlog.logError("SYNTAX ERROR: Unknown identifier: " + token[1] );
    }
    if(interpretAtomReference(token[2], ana, roa) != 0) 
    {
        errlog.logError("SYNTAX ERROR: Unknown identifier: " + token[1] );
    }
    if(interpretAtomReference(token[3], anb, rob) != 0) 
    {
        errlog.logError("SYNTAX ERROR: Unknown identifier: " + token[1] );
    }
    if(interpretAtomReference(token[4], anj, roj) != 0) 
    {
        errlog.logError("SYNTAX ERROR: Unknown identifier: " + token[1] );
    }

    return 0;
}

int DihedralDefinition::printDefinitionLine(FILE * file)
{
    fprintf(file, " IMPROPER ");
    if(roi < 0)
        fprintf(file, "-");
    if(roi > 0)
        fprintf(file, "+");
    fprintf(file, "%s ", ani.c_str());
    if(roa < 0)
        fprintf(file, "-");
    if(roa > 0)
        fprintf(file, "+");
    fprintf(file, "%s ", ana.c_str());
    if(rob < 0)
        fprintf(file, "-");
    if(rob > 0)
        fprintf(file, "+");
    fprintf(file, "%s ", anb.c_str());
    if(roj < 0)
        fprintf(file, "-");
    if(roj > 0)
        fprintf(file, "+");
    fprintf(file, "%s ", anj.c_str());
    fprintf(file, "\n");
    return 0;
}

int MoleculeDefinition::findImproperDef(const DihedralDefinition &query) const
{
    for(unsigned i=0;i<improper.size();i++)
    {
        if(cmpstring(improper[i].ana,query.ana)&&
            cmpstring(improper[i].ani,query.ani)&&
            cmpstring(improper[i].anj,query.anj)&&
            cmpstring(improper[i].anb,query.anb) ) 
            return i;
    }
    return -1;
}

void hackAtomFlags( AtomParameter& atom )
{
    // This is all one giant HACK, all this info needs to be stored in the .ff file...

    atom.setHydrogen(false);
    atom.setCAlpha(false);
    atom.setBackbone(false);
    atom.setCoreBackbone(false);

    if( atom.name[0] == 'H' ) // OOOO, nasty :-D
        atom.setHydrogen(true);

    if( 0 == atom.pdbname.compare("N") )
    {
        atom.setBackbone(true);
        atom.setCoreBackbone(true);
    }
    else if( 0 == atom.pdbname.compare("C") )
    {
        atom.setBackbone(true);
        atom.setCoreBackbone(true);
    }
    else if( 0 == atom.pdbname.compare("CA") )
    {
        atom.setBackbone(true);
        atom.setCoreBackbone(true);
        atom.setCAlpha(true);
    }
    else if( 0 == atom.pdbname.compare("O") )
    {
        atom.setBackbone(true);
        atom.setCoreBackbone(true);
    }
    else if( 0 == atom.pdbname.compare("OXT") )
    {
        atom.setBackbone(true);
        atom.setCoreBackbone(true);
    }
    else if( 0 == atom.pdbname.compare("H") )
    {
        atom.setBackbone(true);
    }
    else if( 0 == atom.pdbname.compare("HA") )
    {
        atom.setBackbone(true);
    }
    else if( 0 == atom.pdbname.compare("1H") )
    {
        atom.setBackbone(true);
    }
    else if( 0 == atom.pdbname.compare("2H") )
    {
        atom.setBackbone(true);
    }
    else if( 0 == atom.pdbname.compare("3H") )
    {
        atom.setBackbone(true);
    }
}

int MoleculeDefinition::readMoleculeBlock(FILE * file, ParseErrorLogger &errlog, FFParamSet &ffps){
    double x, y, z;
    unsigned atnum;
    int typenumber;
    std::string linestring;
    int errorstatus = 0;

    // read line for line
    while(!feof(file)) 
    { 
        if(getSingleLine(file, linestring)!=0) break;
        removecomments(linestring,commentchar);
        errlog.incLineNumber();
        std::vector<std::string> token;
        token = chopstr(linestring,wspace);
        if(token.size() <= 0) continue; // empty line
        std::string command = token[0];

        if(cmpstring(command,"ENDAMINOACID")) 
            break;
        if(cmpstring(command,"ENDMOLECULE")) 
            break;

        //printf("---> %s \n",linestring.c_str());
        if(cmpstring(command,"LOAD")) 
        {
            if(token.size() < 2) 
            {
                errlog.logError("SYNTAX ERROR: Insufficient number of parameters (" + int2str(token.size()) + ")");
                continue;
            }
            std::string MolName;
            MolName = token[1];

            int imol = ffps.findMoleculeType(MolName.c_str());
            if(imol < 0)
            {
                errlog.logError("SYNTAX ERROR: Cannot find molecule name '" + MolName + "'");
                continue;
            }

            (*this) = ffps.molecule[imol];
        } 
        else if(cmpstring(command,"ATOM_APPEND_MODE")) 
        {			
            if(token.size() < 2) 
            {
                errlog.logError("SYNTAX ERROR: 'START' or 'END' identifier expected after 'ATOM_APPEND_MODE' ");
                continue;
            }
            else if( token[1].compare("START") == 0 )
            {
                m_AtomAppendMode = false;
            }
            else if( token[1].compare("END") == 0 )
            {
                m_AtomAppendMode = true;
            }
            else
            {
                errlog.logError("SYNTAX ERROR: 'START' or 'END' identifier expected after 'ATOM_APPEND_MODE' ");
                continue;
            }

        } 
        else if(cmpstring(command,"ATOM")) 
        {
            if(token.size() < 7) 
            {
                errlog.logError("SYNTAX ERROR: Insufficient number of parameters (" + int2str(token.size()) + ")" );
                continue;
            }

            std::string rawname;
            std::string pdbname;
            std::string newname;
            std::string restraint;
            double x,y,z;

            if(!token2name(token[1],rawname))
            {
                errlog.logError("SYNTAX ERROR: Atom rawname is not a valid token");
                errorstatus = -1;
                continue;
            }
            if(!token2name(token[2],pdbname))
            {
                errlog.logError("SYNTAX ERROR: Atom pdbname is not a valid token");
                errorstatus = -1;
                continue;
            }
            if(str2double(token[3],x)!=0)
            {
                errlog.logError("SYNTAX ERROR: Float Expected for x position");
                errorstatus = -1;
                continue;
            }
            if(str2double(token[4],y)!=0)
            {
                errlog.logError("SYNTAX ERROR: Float Expected for y position");
                errorstatus = -1;
                continue;
            }
            if(str2double(token[5],z)!=0)
            {
                errlog.logError("SYNTAX ERROR: Float Expected for z position");
                errorstatus = -1;
                continue;
            }
            newname = token[6]; // 'newname' is allowed characters like '*'

            if(token.size() >= 8)
            {
                restraint = token[7];
            }

            typenumber = ffps.findAtomType(newname);

            if(typenumber < 0) 
            { 
                // newname was not found
                errlog.logError("SYNTAX ERROR: Identifier '" + std::string(&newname[0]) + "' unknown");
                continue;
            }

            // Now create/assign Type to atom parameter

            // The atom paraemter inherits all of it's Type's parameters
            AtomParameter newatomdefinition(ffps.AtomType[typenumber]);
            // plus has individual ones
            newatomdefinition.FFType = typenumber;
            //printf("%d\n",newatomdefinition.Type);

            newatomdefinition.type_name = newname;
            newatomdefinition.i = (char)atom.size(); // Load in all the atom parameters
            newatomdefinition.rawname = rawname;
            newatomdefinition.pdbname = pdbname;
            newatomdefinition.restraint = restraint;

            newatomdefinition.pos().setTo(x,y,z); // Main state - could potentially be changed by the user
            newatomdefinition.posRef().setTo(x,y,z); // The reference state ... inherited property, irrelevent here (user customisable?)
            newatomdefinition.posGeom().setTo(x,y,z); // Geometry state - **MUST** remain read-only / unchanged

            newatomdefinition.parent = this;
            newatomdefinition.parentname = this->name;
            newatomdefinition.parentletter = this->letter;

            hackAtomFlags( newatomdefinition ); // HACK!!

            if( m_AtomAppendMode )
            {
                atom.push_back(newatomdefinition);
            }
            else
            {
                atom.insert(atom.begin(),newatomdefinition);
            }
        } 
        else if(cmpstring(command,"DELATOM")) 
        {
            if(token.size() < 2)
            {
                errlog.logError("SYNTAX ERROR: Insufficient number of parameters (" + int2str( token.size() ) + ")");
                continue;
            }

            std::string atname;
            if(!token2name(token[1],atname))
            {
                errlog.logError("SYNTAX ERROR: Atom altname is not a valid token");
                continue;
            }

            int Atom_i = findAtomRaw(atname);
            if(Atom_i<0)
            {
                errlog.logError("SYNTAX ERROR: Cannot find atom name '" + atname + "' " );
                continue;
            }

            std::vector< AtomParameter >::iterator it = atom.begin();
            it += Atom_i;
            atom.erase(it);
        } 
        else if(cmpstring(command,"DELCONNECT")) 
        {
            if(token.size() < 2)
            {
                errlog.logError("SYNTAX ERROR: Insufficient number of parameters for CONNECT(" + int2str(token.size() ) + "' ");
                errorstatus = -1;
                continue;
            }
            else if(token[1].c_str()[0] == '+')
            {
                frwdName = "";
            }
            else if(token[1].c_str()[0] == '-')
            {
                backName = "";
            }
            else
            {
                errlog.logError("SYNTAX ERROR: CONNECT must be followed by +AtomName or -AtomName");
                errorstatus = -1;
                continue;
            }
        } 
        else if(cmpstring(command,"CONNECT")) 
        {
            if(token.size() < 5)
            {
                errlog.logError("SYNTAX ERROR: Insufficient number of parameters for CONNECT (" + int2str( token.size() ) + ") ");
                errorstatus = -1;
                continue;
            }
            else if(str2double(token[2],x)!=0)
            {
                errlog.logError("SYNTAX ERROR: Float Expected for x position");
                errorstatus = -1;
                continue;
            }
            else if(str2double(token[3],y)!=0)
            {
                errlog.logError("SYNTAX ERROR: Float Expected for y position");
                errorstatus = -1;
                continue;
            }
            else if(str2double(token[4],z)!=0)
            {
                errlog.logError("SYNTAX ERROR: Float Expected for z position");
                errorstatus = -1;
                continue;
            }
            else if(token[1].c_str()[0] == '+')
            {
                frwdName = token[1].substr(1); // skip leading +/-
                frwdpos.setTo(x,y,z);
            }
            else if(token[1].c_str()[0] == '-')
            {
                backName = token[1].substr(1); // skip leading +/-
                backpos.setTo(x,y,z);
            }
            else
            {
                errlog.logError("SYNTAX ERROR: CONNECT must be followed by +AtomName or -AtomName");
                errorstatus = -1;
                continue;
            }

        } 
        else if(cmpstring(command,"IMPROPER")) 
        {
            DihedralDefinition newimpropertype;
            if(newimpropertype.readDefinitionLine(linestring, errlog) == 0) {
                improper.push_back(newimpropertype);
            }
        } 
        else if(cmpstring(command,"DELIMPROPER")) 
        {
            DihedralDefinition newimpropertype;
            if(newimpropertype.readDefinitionLine(linestring, errlog) == 0) {
                int iimp = findImproperDef(newimpropertype);
                if(iimp < 0 ){
                    errlog.logError("SYNTAX ERROR: Cannot find improper definition : " + linestring );
                    continue;
                }
                std::vector< DihedralDefinition >::iterator it = improper.begin();
                it+=iimp;
                improper.erase(it);
            }

        } 
        else if(cmpstring(command,"TORSION")) 
        {
            DihedralDefinition newtorsiontype;
            if(newtorsiontype.readDefinitionLine(linestring, errlog) == 0) {
                torsion.push_back(newtorsiontype);
            }

        } 
        else if(cmpstring(command,"CHARGE")) 
        {
            if(token.size() < 3) {
                errlog.logError("SYNTAX ERROR: CHARGE parse: insufficient number of parameters");
                continue;
            }
            typenumber = findAtomRaw(token[1]);

            if(typenumber < 0) {
                errlog.logError("SYNTAX ERROR: CHARGE parse: Atom name '" + token[1] + "' unknown" );
                continue;
            }

            double charge;
            if(str2double(token[2],charge)!=0){
                errlog.logError("SYNTAX ERROR: CHARGE parse: Float expected for 'charge', found token: '" + token[1] + "' " );
                continue;
            }
            atom[typenumber].charge = charge; // assign charge
        } 
        else if(cmpstring(command,"VDWRADIUS")) 
        {
            if(token.size() < 3) {
                errlog.logError("SYNTAX ERROR: VDWRADIUS parse: Insufficient number of parameters" );
                continue;
            }
            typenumber = findAtomRaw(token[1]);

            if(typenumber < 0) {
                errlog.logError("SYNTAX ERROR: VDWRADIUS parse: Atom name '" + token[1] + "' ");
                continue;
            }

            double radius;
            if(str2double(token[2],radius)!=0){
                errlog.logError("SYNTAX ERROR: VDWRADIUS parse: Float expected for 'radius', found token: '" + token[1] + "' " );
                continue;
            }
            atom[typenumber].radius = radius; // assign radius
        } 
        else if(cmpstring(command,"EPSILON")) 
        {
            if(token.size() < 3) {
                errlog.logError("SYNTAX ERROR: EPSILON parse: Insufficient number of parameters ");
                errorstatus = -1;
                continue;
            }
            typenumber = findAtomRaw(token[1]);

            if(typenumber < 0) {
                errlog.logError("SYNTAX ERROR: EPSILON parse: Atom name '" + token[1] + "' ");
                errorstatus = -1;
                continue;
            }

            double epsilon;
            if(str2double(token[2],epsilon)!=0){
                errlog.logError("SYNTAX ERROR: EPSILON parse: Float expected for 'epsilon', found token: '" + token[1] + "' ");
                continue;
            }
            if(epsilon < 0.0 ){
                errlog.logError("SYNTAX ERROR: EPSILON parse: Epsilon must be positive!");
                continue;
            }
            atom[typenumber].epsilon = sqrt(epsilon * PhysicsConst::kcal2J / PhysicsConst::Na); // assign epsilon
        } 
        else if(cmpstring(command,"MASS")) 
        {
            if(token.size() < 3) {
                errlog.logError("SYNTAX ERROR: MASS parse: insufficient number of parameters"); 
                errorstatus = -1;
                continue;
            }
            typenumber = findAtomRaw(token[1]);

            if(typenumber < 0) {
                errlog.logError("SYNTAX ERROR: MASS parse: Atom name '" + token[1] + "' ");
                continue;
            }

            double mass;
            if(str2double(token[2],mass)!=0){
                errlog.logError("SYNTAX ERROR: MASS parse: Float expected for 'mass', found token: '" + token[1] + "' ");
                continue;
            }
            atom[typenumber].mass = mass * PhysicsConst::amu; // assign mass
        } 
        else 
        {
            // if it's not a special keyword then see if it's a userdefined
            // keyword (a custom property)
            // the syntax is then:
            // NAMEOFPROP ATOMTYPENAME DATA
            int iprop = ffps.findCustomProperty(token[0]);
            if(iprop>=0)
            {
                if(token.size() < 2)
                {
                    errlog.logError("SYNTAX ERROR: Atom type name expected after custom property '" + command + "' ");
                    continue;
                }

                typenumber = findAtomRaw(token[1]);

                if(typenumber < 0)
                {
                    errlog.logError("SYNTAX ERROR: Atom name '" + token[1] + "'  unknown");
                    continue;
                }

                if(token.size() < 3)
                {
                    errlog.logError("SYNTAX ERROR: Data expected after custom property 'i" + command + "' and atom type name" );
                    continue;
                }

                atom[typenumber].addProperty(token[0],token[2]);
            }
            else
            {
                errlog.logError("ERROR: Identifier '" + command + "' unknown"); 
                continue;
            }
        }
    }

    double totalcharge = 0;
    for(atnum = 0; atnum < atom.size(); atnum++)
    {
        decodeRestraintList(atom[atnum].restraint, atom[atnum].r_covalent, errlog);
        totalcharge += atom[atnum].charge;
    }

    // check that total charge is integral
    double chargedev1 = fabs(floor(totalcharge) - totalcharge);
    double chargedev2 = fabs(ceil(totalcharge) - totalcharge);
    double chargedev = min(chargedev1, chargedev2);

    if(chargedev > 0.05) 
    {
        errlog.logWarning("FORCEFIELD WARNING: Residue '" +  std::string(&name[0]) + "' has non-integral total charge " + double2str(totalcharge) );
    }

    //printf("Totalcharge: %8.4lf\n",totalcharge);
    // Now do a check of the restraint information looking for errors
    // in the definition file
    if(checkRestraintList( errlog ) != 0)
    {
        errlog.logError("ERROR: Errors found in restraint lists - check definition files ");
    }

    return errorstatus; //ok , error free?
}

int MoleculeDefinition::writeMoleculeBlock(FILE *file, FFParamSet &ffps)
{
    size_t j,k;
    fprintf(file, "MOLECULE %c %s\n", letter, &name[0]);
    for(j = 0; j < atom.size(); j++)
    {
        if(atom[j].valid == 0)
            continue; // ignore invalid atoms

        fprintf(file, " ATOM %-4s %-4s %6.3lf %6.3lf %6.3lf %-10s ",
            atom[j].rawname.c_str(),
            atom[j].pdbname.c_str(),
            atom[j].pos().x,
            atom[j].pos().y, 
            atom[j].pos().z, 
            &ffps.AtomType[atom[j].FFType].name[0]);

        for(k = 0; k < atom[j].r_covalent.size(); k++)
        {
            if( atom[j].r_covalent[k].i >= 0 )
            {
                if(atom[atom[j].r_covalent[k].i].valid == 0)
                    continue; // ignore invalid atoms
            }

            if(k != 0)
                fprintf(file, ",");
            if(atom[j].r_covalent[k].roi < 0)
                fprintf(file, "-");
            if(atom[j].r_covalent[k].roi > 0)
                fprintf(file, "+");

            fprintf(file, "%s", atom[j].r_covalent[k].ani.c_str());
        }
        fprintf(file, "\n");
    }

    if(0 != frwdName.size())
    {
        fprintf(file, " CONNECT +%s %8.4lf %8.4lf %8.4lf \n",
            frwdName.c_str(),
            frwdpos.x,
            frwdpos.y,
            frwdpos.z);
    }

    if(0 != backName.size())
    {
        fprintf(file, " CONNECT -%s %8.4lf %8.4lf %8.4lf \n",
            backName.c_str(),
            backpos.x,
            backpos.y,
            backpos.z);
    }

    for(j = 0; j < atom.size(); j++)
    {
        // if charge of individual atom differs from atom Type, print explicit charge
        if(atom[j].valid == 0)
            continue; // ignore invalid atoms
        if(atom[j].charge != ffps.AtomType[atom[j].FFType].charge)
        {
            if(atom[j].comment.c_str()=="")
            {
                fprintf(file, " CHARGE %-4s %8.5lf \n",
                    atom[j].rawname.c_str(), atom[j].charge);
            }
            else
            {
                fprintf(file, " CHARGE %-4s %8.5lf # %s\n",
                    atom[j].rawname.c_str(), atom[j].charge, atom[j].comment.c_str());
            }
        }
    }

    for(j = 0; j < atom.size(); j++)
    {
        // if charge of individual atom differs from atom Type, print explicit mass
        if(atom[j].valid == 0)
            continue; // ignore invalid atoms
        if(atom[j].mass != ffps.AtomType[atom[j].FFType].mass)
        {
            fprintf(file, " MASS %-4s %8.5lf \n",
                atom[j].rawname.c_str(), 
                atom[j].mass / PhysicsConst::amu
                );
        }
    }

    for(j = 0; j < atom.size(); j++)
    {
        // if charge of individual atom differs from atom Type, print explicit radius
        if(atom[j].valid == 0)
            continue; // ignore invalid atoms
        if(atom[j].radius != ffps.AtomType[atom[j].FFType].radius)
        {
            fprintf(file, " VDWRADIUS %-4s %8.5lf \n",
                atom[j].rawname.c_str(), atom[j].radius);
        }
    }

    for(j = 0; j < improper.size(); j++)
    {
        improper[j].printDefinitionLine(file);
    }

    fprintf(file, "ENDMOLECULE\n\n\n");

    return 0;
}

bool MoleculeDefinition::hasBackLink() const 
{ 
    return backName.size() != 0; 
}

bool MoleculeDefinition::hasFrwdLink() const 
{ 
    return frwdName.size() != 0; 
}

int Section::readSectionBlock(FILE * file, ParseErrorLogger &errlog)
{
    std::string linestring;
    int errorstatus = 0;

    while(!feof(file))
    { 
        // read line for line
        if(getSingleLine(file, linestring)!=0) break;
        sectionline.push_back(linestring);
        removecomments(linestring,commentchar);

        errlog.incLineNumber();	
        std::vector<std::string> token;
        token = chopstr(linestring,wspace);
        if(token.size() <= 0) 
            continue; // empty line

        std::string command = token[0];
        if(cmpstring(command,"ENDSECTION"))
        {
            sectionline.pop_back(); // remove that last line
            break;
        }
    }

    return errorstatus; //ok , error free
}


/// --------------------------------------------------------------
// FFParamSet Class Function Definitions

void FFParamSet::settodefault()
{
    identgiven = false;
    autobonds = true;
    autoangles = true;
    autotorsions = true;

    ffidentifier = "";

    Vdw14Scaling = 1.0;
    Elec14Scaling = 1.0;

    // Please remember to add to this list any new keywords added to the ff format!
    m_ReservedKeywords.push_back("ALIAS");
    m_ReservedKeywords.push_back("ANGLE");
    m_ReservedKeywords.push_back("ATOM");
    m_ReservedKeywords.push_back("ATOM_APPEND_MODE");
    m_ReservedKeywords.push_back("AUTOANGLES");
    m_ReservedKeywords.push_back("AUTOBONDS");
    m_ReservedKeywords.push_back("AUTOTORSIONS");
    m_ReservedKeywords.push_back("BOND");
    m_ReservedKeywords.push_back("CHARGE");
    m_ReservedKeywords.push_back("CLASS");
    m_ReservedKeywords.push_back("CLASS_CAP");
    m_ReservedKeywords.push_back("CLASS_POSTFIX");
    m_ReservedKeywords.push_back("CLASS_PREFIX");
    m_ReservedKeywords.push_back("CONNECT");
    m_ReservedKeywords.push_back("DELATOM");
    m_ReservedKeywords.push_back("DELCONNECT");
    m_ReservedKeywords.push_back("DELIMPROPER");
    m_ReservedKeywords.push_back("DELTORSION");
    m_ReservedKeywords.push_back("ELEC14SCALING");
    m_ReservedKeywords.push_back("ENDAMINOACID");
    m_ReservedKeywords.push_back("ENDMOLECULE");
    m_ReservedKeywords.push_back("ENDSECTION");
    m_ReservedKeywords.push_back("IDENT");
    m_ReservedKeywords.push_back("IMPROPER");
    m_ReservedKeywords.push_back("INCLUDE");
    m_ReservedKeywords.push_back("INFO");
    m_ReservedKeywords.push_back("LOAD");
    m_ReservedKeywords.push_back("MASS");
    m_ReservedKeywords.push_back("MOLECULE");
    m_ReservedKeywords.push_back("CUSTOMPROPERTY");
    m_ReservedKeywords.push_back("NOAUTOANGLES");
    m_ReservedKeywords.push_back("NOAUTOBONDS");
    m_ReservedKeywords.push_back("NOAUTOTORSIONS");
    m_ReservedKeywords.push_back("SECTION");
    m_ReservedKeywords.push_back("TORSION");
    m_ReservedKeywords.push_back("TYPE");
    m_ReservedKeywords.push_back("VDW14SCALING");
    m_ReservedKeywords.push_back("VDWRADIUS");
}

void FFParamSet::printAtomParameters(AtomParameter * apar)
{
    printf(" Z:%2d %4.1lfamu %4.1lfA(rad) %5.2lf % 5.2lf ",
        apar->Z, apar->mass / PhysicsConst::amu, apar->radius, sqr(apar->epsilon) * PhysicsConst::J2kcal * PhysicsConst::Na, apar->charge);
}

bool FFParamSet::lookupLongName( char _SourceID, std::string& _DestID ) const
{
    int index = findMoleculeType(std::string(1,_SourceID));
    if( index >= 0 )
    {
        _DestID = molecule[index].s_name();
        return true;
    }
    return false;
}

bool FFParamSet::lookupShortName( const std::string& _SourceID, char& _DestID ) const
{
    int index = findMoleculeType(_SourceID);
    if( index >= 0 )
    {
        _DestID = molecule[index].letter;
        return true;
    }
    return false;
}

bool FFParamSet::lookupAlias( const std::string& _SourceID, std::string& _DestID ) const
{	
    int index = findMoleculeType(_SourceID);
    if( index >= 0 )
    {
        _DestID = molecule[index].s_name();
        return true;
    }
    return false;
}

int FFParamSet::findAtomType(const std::string &name) const
{
    unsigned i;
    if(cmpstring(name, "X")) //Found it
        return -1; //-1 means any atom

    for(i = 0; i < AtomType.size(); i++) {
        if(cmpstring(name, AtomType[i].name)) //Found it
            return i; //break out of for loop & return
    }
    return -2; //-2 is an error signal
}

int AtomTypeParameter::AtomTypeReadLine(const std::string linestring, ParseErrorLogger &errlog)
{
    std::vector<std::string> token;
    token = chopstr(linestring,wspace);
    if(token.size() <= 0) return 0; // empty line

    // Warn if parameter list is incomplete
    if(token.size() < 7)
    {
        errlog.logError("SYNTAX ERROR: Definition incomplete (" + int2str(token.size()) + "/5 parameters given) ");
        return -1;
    }

    name = token[1]; // set name of stom Type
    used = 1; // set used to be true;
    if( str2int(token[2],Z) != 0){          errlog.logError("SYNTAX ERROR: Integer expected as parameter 2");return -1; }
    if( str2double(token[3],mass) != 0){    errlog.logError("SYNTAX ERROR: Float expected as parameter 3 "); return -1; }
    if( str2double(token[4],radius) != 0){  errlog.logError("SYNTAX ERROR: Float expected as parameter 4 "); return -1; }
    if( str2double(token[5],epsilon) != 0){ errlog.logError("SYNTAX ERROR: Float expected as parameter 5 "); return -1; }
    if( str2double(token[6],charge) != 0){  errlog.logError("SYNTAX ERROR: Float expected as parameter 6 "); return -1; }

    // Units and other calculations go here !
    // Epsilon is given in kcal/mol; convert it into units of sqrt(J)
    if(epsilon < 0.0 )
    {
        errlog.logError("SYNTAX ERROR: EPSILON parse: Epsilon must be positive! ");
        return -1;
    }
    epsilon = sqrt(epsilon * PhysicsConst::kcal2J / PhysicsConst::Na);

    //The mass is given in atomic mass units (amu) not kg
    mass *= PhysicsConst::amu; // convert to kg !

    return 0;
}

int FFParamSet::findBondType(int ti, int tj) const
{
    int curcand = -1;
    int curwildcards = 4;
    int result;
    unsigned nbond;

    //printf("Looking for %3d %3d Type bond ... \n",ti,tj);

    for(nbond = 0; nbond < BondType.size(); nbond++)
    {
        result = BondType[nbond].cmpBond(ti, tj);
        if(result < 0)
            continue; // no match, go to next

        /* printf(" Found %3d %3d, with %d wildcards\n",
        BondType[nbond].i,
        BondType[nbond].j,result);
        */
        if(result == curwildcards)
        {
            printf("WARNING: Bond Type ambiguity due to wildcards. Ignoring second Type\n");
            continue;
        }
        if(result < curwildcards)
        {
            curwildcards = result;
            curcand = nbond;
        }
    }
    if(curcand < 0)
    {
        printf("WARNING: No appropriate Bond Type (%s %s) found - ignoring bond\n",
            AtomType[ti].name.c_str(), AtomType[tj].name.c_str());
        return -1;
    }

    return curcand;
}



int FFParamSet::findAngleType(int ti, int ta, int tj) const
{
    int curcand = -1;
    int curwildcards = 4;
    int result;
    unsigned nang;

    // printf("Looking for %3d %3d %3d Type angle ... \n",ti,ta,tj);

    for(nang = 0; nang < AngleType.size(); nang++)
    {
        result = AngleType[nang].cmpAngle(ti, ta, tj);
        if(result < 0)
            continue; // no match, go to next

        /* printf(" Found %3d %3d %3d, with %d wildcards\n",
        AngleType[nang].i,
        AngleType[nang].a,
        AngleType[nang].j,result);
        */
        if(result == curwildcards)
        {
            printf(
                "WARNING: Angle Type ambiguity due to wildcards or angle Type duplication\n Looking for angletype %s %s %s Ignoring second Type\n",
                AtomType[ti].name.c_str(), AtomType[ta].name.c_str(), AtomType[tj].name.c_str());
            continue;
        }
        if(result < curwildcards)
        {
            curwildcards = result;
            curcand = nang;
        }
    }
    if(curcand < 0)
    {
        printf("WARNING: No appropriate Angle Type (%s %s %s) found - ignoring angle\n",
            AtomType[ti].name.c_str(), AtomType[ta].name.c_str(), AtomType[tj].name.c_str());
        return -1;
    }

    return curcand;
}


int FFParamSet::findTorsionType(int ti, int ta, int tb, int tj) const
{
    int curcand = -1;
    int curwildcards = 4;
    int result;
    unsigned ntor;

    // printf("Looking for %3d %3d %3d %3d Type torsion ... \n",ti,ta,tb,tj);

    for(ntor = 0; ntor < TorsionType.size(); ntor++)
    {
        result = TorsionType[ntor].cmpTorsion(ti, ta, tb, tj, 0); //0 for torsion as opposed to imnproper (1)
        if(result < 0)
            continue; // no match, go to next

        if(result == curwildcards)
        {
            printf("WARNING: Torsion Type ambiguity due to wildcards. Ignoring second Type\n");
            printf("WARNING: Torsion: %s %s %s %s - wildcards: %d torsion %d %d %d %d\n",
                AtomType[ti].name.c_str(),
                AtomType[ta].name.c_str(),
                AtomType[tb].name.c_str(),
                AtomType[tj].name.c_str(),
                curwildcards,
                TorsionType[ntor].i,
                TorsionType[ntor].a,
                TorsionType[ntor].b,
                TorsionType[ntor].j);
            continue;
        }
        if(result < curwildcards)
        {
            curwildcards = result;
            curcand = ntor;
        }
    }
    if(curcand < 0)
    {
        printf("WARNING: No appropriate Torsion Type (%s %s %s %s)found - ignoring torsion\n",
            AtomType[ti].name.c_str(),
            AtomType[ta].name.c_str(),
            AtomType[tb].name.c_str(),
            AtomType[tj].name.c_str());
        return -1;
    }

    return curcand;
}

// This function looks for the correct improper parameter(s) for a given 4 atom
// types - it uses the same array of classes as findTorsionType because impropers
// and torsions are treated together at this level, but hands a different parameter to
// cmp Torsion (namely 1 instead of 0) to indicate that it's looking for an improper only
int FFParamSet::findImproperType(int ti, int ta, int tb, int tj) const
{
    int curcand = -1;
    int curwildcards = 4;
    int result;
    unsigned ntor;

    // printf("Looking for %3d %3d %3d %3d Type torsion ... \n",ti,ta,tb,tj);

    for(ntor = 0; ntor < TorsionType.size(); ntor++)
    {
        result = TorsionType[ntor].cmpTorsion(ti, ta, tb, tj, 1); //1 for improper as opposed to torsion (1)
        if(result < 0)
            continue; // no match, go to next

        /* printf(" Found %3d %3d %3d %3d, with %d wildcards\n",
        TorsionType[ntor].i,
        TorsionType[ntor].a,
        TorsionType[ntor].b,
        TorsionType[ntor].j,result);
        */
        if(result == curwildcards)
        {
            printf("WARNING: Improper Type ambiguity due to wildcards. Ignoring second Type\n");
            continue;
        }
        if(result < curwildcards)
        {
            curwildcards = result;
            curcand = ntor;
        }
    }
    if(curcand < 0)
    {
        printf("WARNING: No appropriate Improper Type (%s %s %s %s)found - ignoring improper\n",
            AtomType[ti].name.c_str(),
            AtomType[ta].name.c_str(),
            AtomType[tb].name.c_str(),
            AtomType[tj].name.c_str());
        return -1;
    }

    return curcand;
}

int FFParamSet::findMoleculeType_withoutAlias(const std::string &queryname) const
{
    for(unsigned i = 0; i < molecule.size(); i++)
    {
        if(strcmp(&molecule[i].name[0], queryname.c_str()) == 0)
        {
            return i;
        }
    }
    return -1;
}

int FFParamSet::findAliasName(const std::string &_queryAlias,std::string &_returnAliasName) const
{
    for(unsigned i = 0; i < AliasDef.size(); i++)
    {
        if(AliasDef[i].cmpAlias(_queryAlias,_returnAliasName))
        {
            return 0;
        }
    }
    _returnAliasName = std::string();
    return -1;
}

//finds a molecule from the database
int FFParamSet::findMoleculeType(const std::string &MolName) const
{
    int result;

    result = findMoleculeType_withoutAlias(MolName);
    if(result >= 0) return result; // found it without the alias ?

    std::string foundAliasName;
    if(findAliasName(MolName, foundAliasName) < 0) 
    {
        return -1;
    }

    result = findMoleculeType_withoutAlias(foundAliasName);
    if(result >= 0) return result; // found with the alias ?

    return -1; // didn't find any molecule under that name
}

int FFParamSet::findSection(const std::string &queryname) const
{
    for(unsigned i = 0; i < section.size(); i++)
    {
        if(cmpstring(section[i].sectionname, queryname)) return i;
    }
    return -1;
}

int FFParamSet::findCustomProperty(const std::string &queryname) const
{
    for(unsigned i = 0; i < m_CustomProperty.size(); i++)
    {
        if(cmpstring(m_CustomProperty[i], queryname)) return i;
    }
    return -1;
}


int FFParamSet::findReservedKeyword(const std::string &queryname) const
{
    for(unsigned i = 0; i < m_ReservedKeywords.size(); i++)
    {
        if(cmpstring(m_ReservedKeywords[i], queryname)) return i;
    }
    return -1;
}


int FFParamSet::getSection(const std::string &queryname, Section &result) const
{
    int sectionno = findSection(queryname);
    if(sectionno<0)return -1;

    result = section[sectionno];
    return 0;
}


int FFParamSet::readLib(const std::string &filename, int level)
{
    FILE *file;
    double value;
    int errorstatus = 0;

    int ignore = 0;
    std::string ignore_release("");

    std::string linestring;

    std::string fullfilename = LibraryPathStore::getSingleton()->findFullFilename(filename);

    file = fopen(fullfilename.c_str(), "r"); // Open the data file
    if(file == NULL) 
    {
        throw(IOException("Forcefield file '" + filename + "' not found"));
    }

    if( level == 0 )
    {
        printf("Reading library: '%s' \n",fullfilename.c_str());
    }
    else
    {
        printf("  Reading include library: File:'%s' Level:'%d'\n",fullfilename.c_str(),level);
    }

    // creat the error logger

    ParseErrorLogger errlog( fullfilename );
    //Start reading parameter file...

    while(!feof(file))
    {
        // read line for line
        if(getSingleLine(file, linestring)!=0) break;
        removecomments(linestring,commentchar);
        errlog.incLineNumber();
        std::vector<std::string> token;
        token = chopstr(linestring,wspace);
        if(token.size() <= 0) continue; // empty line
        std::string command = token[0];

        //printf("Line %3d: %s\n",line,linestring.c_str());

        if(ignore)
        {
            if(cmpstring(command,ignore_release))
                ignore = false; //end ignore section
            continue;
        }

        //printf("--> %s \n",linestring.c_str());

        if(cmpstring(command,"INFO"))
        {
            // INFO - print trailing text
            printf("%s\n", linestring.c_str() );
        }
        else if(cmpstring(command,"IDENT"))
        {
            // IDENT Forcefield identifier string
            if(token.size() < 2)
            {
                errlog.logError("SYNTAX ERROR: Identifier expected after 'IDENT' ");
                continue;
            }

            ffidentifier = token[1];
            identgiven = true;

        }
        else if(cmpstring(command,"INCLUDE"))
        {
            // INCLUDE a second file
            if(token.size() < 2)
            {
                errlog.logError("SYNTAX ERROR: Identifier expected after 'INCLUDE' ");
                continue;
            }

            StringBuilder sb(filename);
            if( token[1].length() > 2 && token[1][1] == ':' && isalpha(token[1][0]) ) // ':' as in e:\mypath\myfile.ext
            {
                sb.clear();
                sb.setTo(token[1]);
            }
            else
            {
                // make the filename relative to the current name...
                size_t chopper = sb.LastOf("\\/");
                if( chopper == SIZE_T_FAIL )
                {
                    sb.clear();
                }
                else
                {
                    sb.erase(chopper+1,sb.size()-chopper-1);
                }
                sb.append(token[1]);
            }

            // check if file exists
            if(!IO::fileExists(LibraryPathStore::getSingleton()->findFullFilename(sb.toString())))
            {
                errlog.logWarning("WARNING: INCLUDE statement refers to non-existant file '" + sb.toString() + "' - this may or may not be a problem");
                continue;
            }

            // read file
            if( 0 != readLib( LibraryPathStore::getSingleton()->findFullFilename(sb.toString()), level+1 ) )
            {
                errlog.logError("IO ERROR: Could not open included file!");
                continue;
            }
        }
        else if(cmpstring(command,"ALIAS"))
        {
            // ALIAS
            if(token.size() < 3)
            {
                errlog.logError("SYNTAX ERROR: Two identifiers expected after 'ALIAS' ");
                continue;
            }

            std::string aliasSource;
            std::string aliasDest;

            if(!token2name(token[1],aliasSource))
            {
                errlog.logError("SYNTAX ERROR: Residue 'aliasSource'. Syntactically invalid: '" + token[1] );
                continue;
            }
            if(!token2name(token[2],aliasDest))
            {
                errlog.logError("SYNTAX ERROR: Residues 'aliasDest'. Syntactically invalid: '" + token[2] );
                continue;
            }

            if( 0 == aliasSource.compare(aliasDest) )
            {
                errlog.logError("Warning ALIAS token points to itself, ignored! 'i" + aliasSource + "' == '" + aliasDest + "'");
                continue;
            }

            // Allow alias overrides!
            for( int i = ((int)AliasDef.size()-1); i >= 0; i-- )
            {
                if( 0 == AliasDef[i].getAlias().compare(aliasSource) )
                {
                    std::vector<ResidueAliasDefinition>::iterator it = AliasDef.begin();
                    it+=i;
                    AliasDef.erase(it);
                }
            }

            ResidueAliasDefinition newalias;
            newalias.set(aliasSource,aliasDest);
            AliasDef.push_back(newalias);  // save this alias

            // Resolve Aliases to Aliases
            bool changes;
            while( true )
            {
                changes = false;
                for( size_t i = 0; i < AliasDef.size(); i++ )
                {
                    for( size_t j = 0; j < AliasDef.size(); j++ )
                    {
                        if( 0 == AliasDef[i].getName().compare(AliasDef[j].getAlias()) )
                        {
                            // i is pointing to an alias, it needs to instead point to the real name
                            changes = true;
                            AliasDef[i].set(AliasDef[i].getAlias(),AliasDef[j].getName());
                        }
                    }
                }
                if(!changes) break;
            }
        }
        else if(cmpstring(command,"CLASS"))
        {
            // CLASS
            if(token.size() < 3)
            {
                errlog.logError("SYNTAX ERROR: Two identifiers expected after 'CLASS' ");
                continue;
            }
            addClassMember(token[1],token[2]);

        }
        else if(cmpstring(command,"CLASS_CAP"))
        {
            // CLASS_CAP - defines capping residues
            if(token.size() < 3)
            {
                errlog.logError("SYNTAX ERROR: Three identifiers expected after 'CLASS_CAP'  ");
                continue;
            }
            addClassCap(token[1],token[2]);
        }
        else if(cmpstring(command,"CLASS_PREFIX"))
        {
            // CLASS_PREFIX
            if(token.size() < 4)
            {
                errlog.logError("SYNTAX ERROR: Three identifiers expected after 'CLASS_PREFIX' ");
                continue;
            }
            setClassMapperTermini(token[1],
                Library::PrePostFix(token[2],false),
                Library::PrePostFix(token[3],false)
                );

        }
        else if(cmpstring(command,"CLASS_POSTFIX"))
        {
            // CLASS_POSTFIX
            if(token.size() < 4)
            {
                errlog.logError("SYNTAX ERROR: Three identifiers expected after 'CLASS_POSTFIX' ");
                continue;
            }
            setClassMapperTermini(token[1],
                Library::PrePostFix(token[2],true),
                Library::PrePostFix(token[3],true)
                );

        }
        else if(cmpstring(command,"NOAUTOBONDS"))
        {
            // AUTOANGLES - instruct to autogenerate angles
            autobonds = false;
        }
        else if(cmpstring(command,"NOAUTOANGLES"))
        {
            // AUTOANGLES - instruct to autogenerate angles
            autoangles = false;
        }
        else if(cmpstring(command,"NOAUTOTORSIONS"))
        {
            // AUTOTORSIONS - instruct to autogenerate torsions
            autotorsions = false;
        }
        else if(cmpstring(command,"AUTOBONDS"))
        {
            // AUTOANGLES - instruct to autogenerate angles
            autobonds = true;
        }
        else if(cmpstring(command,"AUTOANGLES"))
        {
            // AUTOANGLES - instruct to autogenerate angles
            autoangles = true;
        }
        else if(cmpstring(command,"AUTOTORSIONS"))
        {
            // AUTOTORSIONS - instruct to autogenerate torsions
            autotorsions = true;
        }
        else if(cmpstring(command,"VDW14SCALING"))
        {
            if(token.size() < 2)
            {
                errlog.logError("SYNTAX ERROR: Identifier expected after '" + command + "' " );
                continue;
            }

            if(str2double(token[1],value)!=0)
            {
                errlog.logError("SYNTAX ERROR: Identifier expected after '" + command + "' " );
                continue;
            }
            Vdw14Scaling = value;
        }
        else if(cmpstring(command,"ELEC14SCALING"))
        {
            if(token.size() < 2)
            {
                errlog.logError("SYNTAX ERROR: Identifier expected after '" + command + "' " );
                continue;
            }

            if(str2double(token[1],value)!=0)
            {
                errlog.logError("SYNTAX ERROR: Floating point number expected after '" + command + "' " );
                continue;
            }
            Elec14Scaling = value;
        }        
        else if(cmpstring(command,"CUSTOMPROPERTY")) 
        {
            // Adding a new custom property
            if(token.size() < 2)
            {
                errlog.logError("SYNTAX ERROR: Identifier (The name of the custom property) expected after '" +  command + "'");
                continue;
            }

            int iprop = findReservedKeyword(token[1]);
            if(iprop>=0)
            {
                errlog.logError("SYNTAX ERROR: Custom property '" + command + "' cannot use reserved keyword - please choose a different name");
            }
            else
            {
                iprop = findCustomProperty(token[1]);
                if(iprop>=0)
                {
                    errlog.logError("SYNTAX ERROR: Custom property '" + command + "' cannot use reserved keyword - please choose a different name");
                }
                else
                {
                    m_CustomProperty.push_back(token[1]);
                }
            }
        }
        else if(cmpstring(command,"SECTION")) 
        {
            // ADDITIONAL SECTION syntax reading/checking
            if(token.size() < 2)
            {
                errlog.logError("SYNTAX ERROR: Identifier expected after '" + command + "' " );
                continue;
            }

            int isection = findSection(token[1]);
            if(isection>=0)
            {
                printf("  Adding to Section '%s'\n",token[1].c_str());
                section[isection].readSectionBlock(file, errlog);
            }
            else
            {
                printf("  Creating new section '%s'\n",token[1].c_str());
                Section newsection(token[1],filename,errlog.getLineNumber() );
                newsection.readSectionBlock(file, errlog);
                section.push_back(newsection);
            }


        }
        else if(cmpstring(command,"TYPE"))
        {
            // TYPE syntax reading/checking
            AtomTypeParameter newatomtype;

            if(newatomtype.AtomTypeReadLine(linestring.c_str(), errlog)!=0)
            {
                errlog.logError("ERROR: Error reading Type entry");
            }
            else
            {
                AtomType.push_back(newatomtype);
            }

        }
        else if(cmpstring(command,"BOND"))
        {
            // BOND syntax reading/checking
            BondTypeParameter newbondtype;
            if(newbondtype.readDefinitionLine(linestring.c_str(),  *this, errlog) == 0)
            {
                BondType.push_back(newbondtype);
            }
            else
            {
                errlog.logError("ERROR: Error reading bond type entry");
            }

        }
        else if(cmpstring(command,"ANGLE"))
        {
            // ANGLE syntax reading/checking
            AngleTypeParameter newangletype;
            if(newangletype.readDefinitionLine(linestring.c_str(),  *this, errlog) == 0)
            {
                AngleType.push_back(newangletype);
            }
            else
            {
                errlog.logError("ERROR: Error reading angle type entry");
            }

        }
        else if((cmpstring(command,"TORSION")) || (cmpstring(command,"IMPROPER")))
        {
            // TORSION/IMPROPER syntax reading/checking
            TorsionTypeParameter newtorsiontype;
            int treturn = newtorsiontype.readDefinitionLine(TorsionType, TorsionType.size(), linestring, *this, errlog);
            if(treturn==0)
            {
                TorsionType.push_back(newtorsiontype);
            }
            else if (treturn<0)
            {
                errlog.logError("ERROR: Error reading torsion type entry");
            }
        } 
        else if(cmpstring(command,"MOLECULE")) 
        { 
            // MOLECULE syntax reading/checking

            if(token.size() < 4) 
            {
                errlog.logWarning("SYNTAX WARNING: Three identifiers expected after 'MOLECULE'\n - Token 1 is the single letter name\n - Token 2 is the three letter name\n - Token 3 is the full name");
            }
            else if( token[1].size() != 1 )
            {
                errlog.logWarning("SYNTAX WARNING: 'MOLECULE' Single letter name token much be of length 1 - truncating!!");
            }
            else if( token[2].size() > 3 )
            {
                errlog.logWarning("SYNTAX WARNING: 'MOLECULE' Three letter name token must be less than length 3 - truncating!!!");
            }

            MoleculeDefinition newmoleculetype; // Make ourselves new molecule
            if(token.size() >=4)
            { 
                // NOTE: These MUST be added AFTER readMoleculeBlock(), as it can call 'LOAD' which overwrites existing data
                newmoleculetype.letter = token[1].c_str()[0]; // add single letter
                if(!token2name(token[2],newmoleculetype.l3_name))
                {
                    errlog.logError("SYNTAX ERROR: Residue 3-letter name '" + newmoleculetype.name + "' is not a valid identifier \n");
                    continue;
                }
                if(!token2name(token[3],newmoleculetype.name))
                {
                    errlog.logError("SYNTAX ERROR: Residue full name '" + newmoleculetype.name + "' is not a valid identifier \n");
                    continue;
                }
                if( newmoleculetype.l3_name.size() > 3 )
                {
                    newmoleculetype.l3_name = newmoleculetype.l3_name.substr(0,3);
                }
            }
            else if(token.size() >=3)
            {
                // NOTE: These MUST be added AFTER readMoleculeBlock(), as it can call 'LOAD' which overwrites existing data
                newmoleculetype.letter = token[1].c_str()[0]; // add single letter
                if(!token2name(token[2],newmoleculetype.name))
                {
                    errlog.logError("SYNTAX ERROR: Residue 3-letter name '" + newmoleculetype.name + "' is not a valid identifier \n");
                    continue;
                }
                // create 3 letter name pointer
                int alen;
                alen = (int) strlen(&newmoleculetype.name[0]) - 3;
                if(alen < 0) alen = 0;
                newmoleculetype.l3_name = &newmoleculetype.name[alen];
            } 

            // read atom information (ENDMOLECULE terminates)
            if(newmoleculetype.readMoleculeBlock(file, errlog, *this) != 0) 
            { 				
                errorstatus = -1;
                continue;
            } 

            if(token.size() >=4)
            { 
                // NOTE: These MUST be added AFTER readMoleculeBlock(), as it can call 'LOAD' which overwrites existing data
                newmoleculetype.letter = token[1].c_str()[0]; // add single letter
                if(!token2name(token[2],newmoleculetype.l3_name))
                {
                    errlog.logError("SYNTAX ERROR: Residue 3-letter name '" + newmoleculetype.name + "' is not a valid identifier \n");
                    continue;
                }
                if(!token2name(token[3],newmoleculetype.name))
                {
                    errlog.logError("SYNTAX ERROR: Residue full name '" + newmoleculetype.name + "' is not a valid identifier \n");
                    continue;
                }
                if( newmoleculetype.l3_name.size() > 3 )
                {
                    newmoleculetype.l3_name = newmoleculetype.l3_name.substr(0,3);
                }
            }
            else if(token.size() >=3)
            {
                // NOTE: These MUST be added AFTER readMoleculeBlock(), as it can call 'LOAD' which overwrites existing data
                newmoleculetype.letter = token[1].c_str()[0]; // add single letter
                if(!token2name(token[2],newmoleculetype.name))
                {
                    errlog.logError("SYNTAX ERROR: Residue 3-letter name '" + newmoleculetype.name + "' is not a valid identifier \n");
                    continue;
                }
                // create 3 letter name pointer
                int alen;
                alen = (int) strlen(&newmoleculetype.name[0]) - 3;
                if(alen < 0) alen = 0;
                newmoleculetype.l3_name = &newmoleculetype.name[alen];
            } 
            molecule.push_back(newmoleculetype);
        } 
        else 
        {
            // if it's not a special keyword then see if it's a userdefined
            // keyword (a custom property)
            // the syntax is then:
            // NAMEOFPROP ATOMTYPENAME DATA
            int iprop = findCustomProperty(token[0]);
            if(iprop>=0)
            {
                if(token.size() < 2) 
                {
                    errlog.logError("SYNTAX ERROR: Atom type name expected after custom property '" + command + "' " );
                    continue;
                }
                int itype = findAtomType(token[1]);
                if(itype < 0)
                {
                    errlog.logError("SYNTAX ERROR: Unknown atom type '" + token[1] + "' ");
                    continue;
                }
                if(token.size() < 3) 
                {
                    errlog.logError("SYNTAX ERROR: Data expected after custom property '" + command + "' and atom type name ");
                    continue;
                }

                AtomType[itype].addProperty(token[0],token[2]);
            }
            else
            {
                errlog.logError("SYNTAX ERROR: Unknown identifier '" + command + "' ");
                continue;
            }
        }
    }

    // final checks
    if(ignore) {
        errlog.logError("ERROR: Unexpected end of file - no '" + ignore_release + "' found ");
    }

    printf("  Finished reading Forcefield Definition File %s ", filename.c_str());
    if(  (errorstatus == 0) && ( errlog.getErrorCount() == 0 ))
    {
        printf("no errors\n");
    }
    else
    {
        printf("errors occured\n");
        if(level==0) throw(ParseException("Errors occured in forcefield file(s) - see above for details"));
    }

    // A little code used to make variations of forcefields
    /*

    printf("WARNING: UNIFYING FORCEFIELD !!!! \n");
    unifyForcefield();
    checkUsedTypes();

    // print all gathered information
    char newfilename[256];
    strcpy(&newfilename[0],filename);
    strcat(&newfilename[0],".new");
    write(&newfilename[0]);
    return -1;
    */

    fclose(file);
    return errorstatus;
}


int FFParamSet::unifyForcefield()
{
    unsigned i, j, k;
    int hcount;
    // New Type united atom forcefield

    // add atomtypes CT1,CT2,CT3
    /*
    strcpy(&AtomType[AtomType.size()].name[0], "CT1");
    AtomType[AtomType.size()].p = AtomType[0].p;
    AtomType[AtomType.size()].used = 1;
    AtomType.size()++;
    strcpy(&AtomType[AtomType.size()].name[0], "CT2");
    AtomType[AtomType.size()].p = AtomType[0].p;
    AtomType[AtomType.size()].used = 1;
    AtomType.size()++;
    strcpy(&AtomType[AtomType.size()].name[0], "CT3");
    AtomType[AtomType.size()].p = AtomType[0].p;
    AtomType[AtomType.size()].used = 1;
    AtomType.size()++;
    */
    for(i = 0; i < molecule.size(); i++)
    {
        // go through all molecules
        for(j = 0; j < molecule[i].atom.size(); j++)
        {
            // and their atoms
            //if( molecule[i].atom[j].Z != 6 ) continue;// if not carbon go to next atom
            // find all the hydrogens

            // exclusion rules - HA, aromatics ?
            /*
            if((strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A1") != 0) &&
            (strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A9") != 0) &&
            (strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A42") != 0) &&
            (strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A65") != 0) &&
            (strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A66") != 0) &&
            (strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A67") != 0) &&
            (strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A68") != 0) &&
            (strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A69") != 0) &&
            (strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A70") != 0) &&
            (strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A71") != 0) &&
            (strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A72") != 0) &&
            (strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A73") != 0) &&
            (strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A74") != 0) &&
            (strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A75") != 0) &&
            (strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A76") != 0) &&
            (strcmp(&AtomType[molecule[i].atom[j].FFType].name[0], "A77") != 0) )
            {
            continue;
            }

            */

            //			if(strcmp(molecule[i].atom[j].pdbname.c_str(), "CA") == 0) {
            //				molecule[i].atom[j].radius += 0.00087; // sort of mark them
            //				continue;
            //			}

            hcount = 0;
            for(k = 0; k < molecule[i].atom[j].r_covalent.size(); k++)
            {
                // scan through all the covalently bound atoms
                if(molecule[i].atom[j].r_covalent[k].i < 0)
                    continue; //ignore if linked to atom in other residue

                if(molecule[i].atom[molecule[i].atom[j].r_covalent[k].i].Z == 1)
                {
                    // found a hydrogen
                    // now add charge of the bound hydrogens to the carbon atom
                    molecule[i].atom[j].charge += molecule[i].atom[molecule[i].atom[j].r_covalent[k].i].charge;
                    // also add the masses of the hydrogens to the mother atom
                    molecule[i].atom[j].mass += molecule[i].atom[molecule[i].atom[j].r_covalent[k].i].mass;

                    molecule[i].atom[molecule[i].atom[j].r_covalent[k].i].valid = 0; // invalidate the hydrogen
                    hcount++;
                }
            }
            //molecule[i].atom[j].Type = (char)(AtomType.size() - hcount);
        }
    }
    return 0;
}

//this reduces the forecefield parameter set to the minimum, by analysing whihc
//atom types are actually used in the molecule definitions, and according to that
//cutting out unused atom types and covalent structures (angles/bonds/torsions/impropers)
//that contain them. This is done by marking them as unused (by setting used to 0)
//writeForcefieldParaemterfile(...) then ignores these during print.

void FFParamSet::checkUsedTypes()
{
    unsigned t, i, j;

    //check filetypes first
    for(t = 0; t < AtomType.size(); t++)
    {
        AtomType[t].used = 0;
        // look through all molecule definitions
        for(i = 0; i < molecule.size(); i++)
        {
            // go through all aminoacid types
            for(j = 0; j < molecule[i].atom.size(); j++)
            {
                // an their atoms
                if(molecule[i].atom[j].valid == 0)
                    continue; // ignore invalid atoms
                if(molecule[i].atom[j].FFType == t)
                {
                    AtomType[t].used = 1; // mark as used
                    break;
                }
            }
            if(AtomType[t].used == 1)
                break; // dont bother looking further if found
        }
    }

    for(i = 0; i < BondType.size(); i++)
    {
        if(BondType[i].i != -1)
            if(AtomType[BondType[i].i].used == 0)
                BondType[i].used = 0; // set used flag to 0 if atomtype is unused
        if(BondType[i].j != -1)
            if(AtomType[BondType[i].j].used == 0)
                BondType[i].used = 0; // set used flag to 0 if atomtype is unused
    }

    for(i = 0; i < AngleType.size(); i++)
    {
        if(AngleType[i].i != -1)
            if(AtomType[AngleType[i].i].used == 0)
                AngleType[i].used = 0; // set used flag to 0 if atomtype is unused
        if(AngleType[i].a != -1)
            if(AtomType[AngleType[i].a].used == 0)
                AngleType[i].used = 0; // set used flag to 0 if atomtype is unused
        if(AngleType[i].j != -1)
            if(AtomType[AngleType[i].j].used == 0)
                AngleType[i].used = 0; // set used flag to 0 if atomtype is unused
    }
    for(i = 0; i < TorsionType.size(); i++)
    {
        if(TorsionType[i].i != -1)
            if(AtomType[TorsionType[i].i].used == 0)
                TorsionType[i].used = 0; // set used flag to 0 if atomtype is unused
        if(TorsionType[i].a != -1)
            if(AtomType[TorsionType[i].a].used == 0)
                TorsionType[i].used = 0; // set used flag to 0 if atomtype is unused
        if(TorsionType[i].b != -1)
            if(AtomType[TorsionType[i].b].used == 0)
                TorsionType[i].used = 0; // set used flag to 0 if atomtype is unused
        if(TorsionType[i].j != -1)
            if(AtomType[TorsionType[i].j].used == 0)
                TorsionType[i].used = 0; // set used flag to 0 if atomtype is unused
    }
    return;
}


int FFParamSet::writeLib(const std::string &filename)
{
    FILE *file;
    unsigned i;
    int h;
    char *X = "X";
    char *ci, *ca, *cb, *cj;

    file = fopen(filename.c_str(), "w"); // Open the data file
    if(file == NULL)
    {
        printf("Cannot open forcefield file for writing\n");
        return -1;
    }

    fprintf(file, "# Atom Types\n#\n");
    fprintf(file, "# Groups .---<atom>----. .-----vdw-------. .--elec--. \n");
    fprintf(file, "# Type (Z) AtMass Radius Epsilon Charge \n");

    for(i = 0; i < AtomType.size(); i++)
    {
        if(AtomType[i].used == 0)
            continue; // skip unused atom types
        fprintf(file, "TYPE %-10s%2d %7.2lf %7.4lf %7.4lf %7.4lf\n",
            &AtomType[i].name[0],
            AtomType[i].Z,
            AtomType[i].mass / PhysicsConst::amu,
            AtomType[i].radius, sqr(AtomType[i].epsilon) / (PhysicsConst::kcal2J / PhysicsConst::Na), AtomType[i].charge);
    }

    fprintf(file, "\n\n\n\n");

    fprintf(file, "# Bond types:\n");
    fprintf(file, "# Atom A B length(A) force constant(kcal/mol/A^2)\n");
    fprintf(file, "#\n");

    for(i = 0; i < BondType.size(); i++)
    {
        if(BondType[i].used == 0)
            continue; // skip unused types

        if(BondType[i].i == -1)
            ci = X;
        else
            ci = &AtomType[BondType[i].i].name[0];
        if(BondType[i].j == -1)
            cj = X;
        else
            cj = &AtomType[BondType[i].j].name[0];

        fprintf(file, "BOND %3s %3s %7.4lf %7.2lf  %4.1f \n",
            ci, cj, BondType[i].length, BondType[i].forceconstant * PhysicsConst::J2kcal * PhysicsConst::Na,  BondType[i].bondorder );
    }
    fprintf(file, "\n\n\n\n");

    fprintf(file, "# Angle types:\n");
    fprintf(file, "# Atom A B C angle(deg) force constant(kcal/mol/A^2)\n");
    fprintf(file, "#\n");

    for(i = 0; i < AngleType.size(); i++) {
        if(AngleType[i].used == 0)
            continue; // skip unused types

        if(AngleType[i].i == -1)
            ci = X;
        else
            ci = &AtomType[AngleType[i].i].name[0];
        if(AngleType[i].a == -1)
            ca = X;
        else
            ca = &AtomType[AngleType[i].a].name[0];
        if(AngleType[i].j == -1)
            cj = X;
        else
            cj = &AtomType[AngleType[i].j].name[0];

        fprintf(file, "ANGLE %3s %3s %3s %6.2lf %7.2lf \n",
            ci, ca, cj, AngleType[i].angle * Maths::MathConst::OneEightyOverPI, AngleType[i].forceconstant * PhysicsConst::J2kcal * PhysicsConst::Na);
    }
    fprintf(file, "\n\n\n\n");



    fprintf(file, "# Torsion types:\n");
    fprintf(file, "# Atom A B C D Vm(kcal/mol) gamma(deg) n \n");
    fprintf(file, "#\n");

    for(i = 0; i < TorsionType.size(); i++) {
        for(h = 0; h < TorsionType[i].terms; h++) {
            if(TorsionType[i].used == 0)
                continue; // skip unused types
            if(TorsionType[i].Type != 0)
                continue;
            if(TorsionType[i].i == -1)
                ci = X;
            else
                ci = &AtomType[TorsionType[i].i].name[0];
            if(TorsionType[i].a == -1)
                ca = X;
            else
                ca = &AtomType[TorsionType[i].a].name[0];
            if(TorsionType[i].b == -1)
                cb = X;
            else
                cb = &AtomType[TorsionType[i].b].name[0];
            if(TorsionType[i].j == -1)
                cj = X;
            else
                cj = &AtomType[TorsionType[i].j].name[0];

            fprintf(file, "TORSION %3s %3s %3s %3s %6.3lf %6.1lf %4.2lf\n",
                ci, ca, cb, cj,
                TorsionType[i].Vn[h] * PhysicsConst::J2kcal * PhysicsConst::Na, TorsionType[i].gamma[h] * Maths::MathConst::OneEightyOverPI, TorsionType[i].n[h]);
        }
    }
    fprintf(file, "\n\n\n\n");


    fprintf(file, "# Improper types:\n");
    fprintf(file, "# Atom A B C D Vm(kcal/mol) gamma(deg) n \n");
    fprintf(file, "#\n");

    for(i = 0; i < TorsionType.size(); i++) {
        for(h = 0; h < TorsionType[i].terms; h++) {
            if(TorsionType[i].used == 0)
                continue; // skip unused types
            if(TorsionType[i].Type == 0)
                continue; // skip torsions
            if(TorsionType[i].i == -1)
                ci = X;
            else
                ci = &AtomType[TorsionType[i].i].name[0];
            if(TorsionType[i].a == -1)
                ca = X;
            else
                ca = &AtomType[TorsionType[i].a].name[0];
            if(TorsionType[i].b == -1)
                cb = X;
            else
                cb = &AtomType[TorsionType[i].b].name[0];
            if(TorsionType[i].j == -1)
                cj = X;
            else
                cj = &AtomType[TorsionType[i].j].name[0];

            fprintf(file, "IMPROPER %3s %3s %3s %3s %6.3lf %6.1lf %4.2lf\n",
                ci, ca, cb, cj,
                TorsionType[i].Vn[h] * PhysicsConst::J2kcal * PhysicsConst::Na, TorsionType[i].gamma[h] * Maths::MathConst::OneEightyOverPI, TorsionType[i].n[h]);
        }
    }
    fprintf(file, "\n\n\n\n");


    fprintf(file, "#Now define residues/molecule types explicitly\n");
    fprintf(file, "#\n");
    fprintf(file, "# atnam x y z Type Connectivity\n");

    for(i = 0; i < molecule.size(); i++) {
        molecule[i].writeMoleculeBlock(file, *this);
    }



    fclose(file);

    return 0;
}



int FFParamSet::readCHARMMprm(const std::string &filename)
{
    FILE *file;
    int line;
    // double value; C4101 - Unreferenced local variable
    int errorstatus = 0;
    bool debugmode = false;
    int ignore = 0;
    std::string ignore_release("");

    std::string linestring;
    std::string rawlinestring;
    std::string comment;
    std::string mode = ""; // "bonds", "angles", "torsions", ..

    std::string wrapline = "";
    ParseErrorLogger errlog( filename );

    file = fopen(filename.c_str(), "r"); // Open the data file
    if(file == NULL) 
    {
        printf("Forcefieldfile '%s' not found\n", filename.c_str());
        return -1;
    }

    MoleculeDefinition newmoldef;
    bool loadedresi = false;

    //Start reading parameter file
    line = -1;
    printf("INFO Converted from '%s' by pd CHARMM importer \n",filename.c_str());
    std::string beforeline = "";
    while(!feof(file)) 
    { 
        // read line for line

        if(getSingleLine(file, rawlinestring)!=0) break;
        linestring = beforeline + " " + rawlinestring;
        removecomments(linestring,"\12\15");

        if(debugmode) printf("### %s \n",linestring.c_str());

        comment="";
        removecomments(linestring,"!",comment);
        line++;

        wrapline = "";
        std::vector<std::string> token;
        token = chopstr(linestring,wspace);

        if(token.size() <= 0)
        {
            printf("# %s \n",comment.c_str());
            continue; // empty line
        }

        if(cmpstring(token[token.size()-1],"-"))
        {
            beforeline = linestring;
            continue;
        }
        else
        {
            beforeline = "";
        }

        std::string command = token[0];

        if(ignore)
        {
            if(cmpstring(command,ignore_release))
                ignore = false; //end ignore section
            continue;
        }
        if(cmpstring(command,"*")) 
        { 
            // INFO - print trailing text
            printf("  INFO %s\n", linestring.c_str() );
        }
        // ignore certain tags - no equivalent present in pd
        else if(cmpstring(command,"DECL"))
        {
            // no declarations needed in pd
        }
        else if(cmpstring(command,"BONDS"))
        {
            mode = "BONDS";
        }
        else if(cmpstring(command,"ANGLES")||
            cmpstring(command,"THETAS"))
        {
            mode = "ANGLES";
        }
        else if( cmpstring(command,"DIHEDRALS")||
            cmpstring(command,"PHI") ||
            cmpstring(command,"DIHE") )
        {
            mode = "DIHEDRALS";
        }
        else if(!cmpstring(mode,"RESI")&&(
            cmpstring(command,"IMPROPER")||
            cmpstring(command,"IMPHI"))) 
        {
            mode = "IMPROPER";
        }
        else if(cmpstring(command,"MASS"))
        {
            if(token.size() < 4)
            {
                printf("## SYNTAX ERROR: at least 4 identifiers expected after MASS \n");
                printf("## %s",rawlinestring.c_str());
            }

            if(token.size() < 3) continue;
            int itype = findAtomType(token[2]);

            if(itype<0)
            {
                printf("## WARNING Previously unknown atom Type: '%s' - creating new atom Type \n",token[2].c_str());
                printf("## without nonbonded parameters! \n");

                double lZ=0;
                double lMass=1.000;
                switch( token[2].c_str()[0] ){
                    case 'H': lZ = 1; lMass=1.000; break;
                    case 'C': lZ = 6; lMass=12.000; break;
                    case 'N': lZ = 7; lMass=14.000; break;
                    case 'O': lZ = 8; lMass=16.000; break;
                    case 'S': lZ = 16; lMass=32.000; break;
                    case 'P': lZ = 15; lMass=31.000; break;
                };

                AtomTypeParameter newatomtype;
                newatomtype.name = token[2]; // set name of stom Type
                newatomtype.used = 1; // set used to be true;
                newatomtype.Z = int(lZ);
                newatomtype.mass = lMass;
                newatomtype.epsilon = 0;
                newatomtype.mass = 0;
                AtomType.push_back(newatomtype);
                itype = AtomType.size() - 1;
            }

            int lZ;
            double lMass;

            if(token.size() >= 4)
            {
                if( str2double(token[3],lMass) != 0)
                {
                    printf("## SYNTAX ERROR: float expected for Z \n");
                    printf("## %s",rawlinestring.c_str());
                    continue;
                }
                AtomType[itype].mass = lMass * PhysicsConst::amu; // convert to kg !;
            }
            if(token.size() >= 5)
            {
                lZ = ElementSymbol2Z(token[4]);
                if(lZ < 0)
                {
                    printf("## UNKNOWN ELEMENT: '%s' \n", token[4].c_str() );
                }
                else
                {
                    AtomType[itype].Z = lZ;
                }
            }

            printf("TYPE %-10s%2d %7.2lf %7.4lf %7.4lf %7.4lf ",
                AtomType[itype].name.c_str(),
                AtomType[itype].Z,
                AtomType[itype].mass / PhysicsConst::amu,
                AtomType[itype].radius,
                sqr(AtomType[itype].epsilon) / (PhysicsConst::kcal2J / PhysicsConst::Na),
                AtomType[itype].charge);
        }
        else if(cmpstring(command,"NONBONDED"))
        {
            mode = "NONBONDED";
            printf("### NON BONDED Parameters: %s",linestring.c_str());

            // parse this: ?
            //nbxmod 5 atom cdiel ForceSwitch vatom vdistance vfswitch -cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14ac 1.0 wmin 1.5
            /*size_t it;
            for(it=1;it<token.size();it++)
            {
            }*/
        }
        else if(cmpstring(command,"RESI") || cmpstring(command,"RESIDUE"))
        {
            if(token.size() < 2)
            {
                printf("## SYNTAX ERROR: 1 identifier expected after RESI (at least a name) \n");
                printf("## %s",rawlinestring.c_str());
                continue;
            }
            if(loadedresi)
            {
                newmoldef.decodeRestraintList(errlog);
                newmoldef.writeMoleculeBlock(stdout, *this);
                loadedresi = false;
            }
            loadedresi = true;

            newmoldef = MoleculeDefinition();
            mode = "RESI";
            newmoldef.name = token[1];
        }
        else if(cmpstring(command,"PRES")) 
        { 
            // PATCH - currently not implemented in pd
            if(loadedresi)
            {
                newmoldef.decodeRestraintList(errlog);
                newmoldef.writeMoleculeBlock(stdout, *this);
                loadedresi = false;
            }
            mode = "PRES";
        }
        else
        {
            if(cmpstring(mode,"BONDS"))
            {
                if( token.size() < 4 )
                {
                    printf("## SYNTAX ERROR: 4 identifiers expected in BONDS section \n");
                    printf("## %s",rawlinestring.c_str());
                    continue;
                }
                printf("BOND %8s %8s %8s %8s ", token[0].c_str(), token[1].c_str(), token[3].c_str(), token[2].c_str() );
            }
            else if(cmpstring(mode,"ANGLES"))
            { 
                if( token.size() >= 6 )
                {
                    printf("## WARNING: pd does not currently support Urey - Bradley Angle Terms \n");
                    printf("## %s",rawlinestring.c_str());
                }
                else if( token.size() < 5 )
                {
                    printf("## SYNTAX ERROR: at least 5 identifiers expected in ANGLES section \n");
                    printf("## %s", rawlinestring.c_str());
                    continue;
                }
                else
                {
                    printf("ANGLE %8s %8s %8s %8s %8s ",
                        token[0].c_str(),
                        token[1].c_str(),
                        token[2].c_str(),
                        token[4].c_str(),
                        token[3].c_str()
                        );
                }
            }
            else if(cmpstring(mode,"DIHEDRALS"))
            { 
                if( token.size() < 7 )
                {
                    printf("## SYNTAX ERROR: 7 identifiers expected in DIHEDRALS section \n");
                    printf("## %s",rawlinestring.c_str());
                    continue;
                }
                printf("TORSION %8s %8s %8s %8s %8s %8s %8s ",
                    token[0].c_str(), token[1].c_str(), token[2].c_str(), token[3].c_str(),
                    token[4].c_str(), token[6].c_str(), token[5].c_str());
            }
            else if(cmpstring(mode,"IMPROPER"))
            { 
                if( token.size() < 7 )
                {
                    printf("## SYNTAX ERROR: 7 identifiers expected in IMPROPER section: \n");
                    printf("## %s",rawlinestring.c_str());
                    continue;
                }
                printf("IMPROPER %8s %8s %8s %8s %8s %8s %8s ",
                    token[0].c_str(), token[1].c_str(), token[2].c_str(), token[3].c_str(),
                    token[4].c_str(), token[6].c_str(), token[5].c_str());
            }
            else if(cmpstring(mode,"NONBONDED"))
            { 
                if( token.size() < 4 )
                {
                    printf("## SYNTAX ERROR: at least 4 identifiers expected in NONBONDED section: \n");
                    printf("## %s",rawlinestring.c_str());
                    continue;
                }
                //printf("TYPE %8s ",token[0].c_str());
                // Guess Z and Mass:
                double lZ=0;
                double lMass=1.000;
                //printf(" '%c' \n",token[0].c_str()[0] );
                switch( token[0].c_str()[0] ){
                    case 'H': lZ = 1; lMass=1.000; break;
                    case 'C': lZ = 6; lMass=12.000; break;
                    case 'N': lZ = 7; lMass=14.000; break;
                    case 'O': lZ = 8; lMass=16.000; break;
                    case 'S': lZ = 16; lMass=32.000; break;
                    case 'P': lZ = 15; lMass=31.000; break;
                };

                //printf("%2d %7.4lf ",int(lZ),lMass);
                //printf(" %8s %8s ", token[3].c_str(), token[2].c_str());

                AtomTypeParameter newatomtype;

                newatomtype.name = token[0]; // set name of stom Type
                newatomtype.used = 1; // set used to be true;
                newatomtype.Z = int(lZ);
                newatomtype.mass = lMass;

                if( str2double(token[3],newatomtype.radius) != 0)
                {
                    printf("## SYNTAX ERROR: float expected for radius \n");
                    printf("## %s",rawlinestring.c_str());
                    continue;
                }
                if( str2double(token[2],newatomtype.epsilon) != 0)
                {
                    printf("## SYNTAX ERROR: float expected for epsilon \n");
                    printf("## %s",rawlinestring.c_str());
                    continue;
                }

                // Units and other calculations go here !
                //epsilon is given in kcal/mol; convert it into units of sqrt(J)
                newatomtype.epsilon = sqrt(-newatomtype.epsilon * PhysicsConst::kcal2J / PhysicsConst::Na);

                //The mass is given in atomic mass units (amu) not kg
                newatomtype.mass *= PhysicsConst::amu; // convert to kg !

                AtomType.push_back(newatomtype);
            } 
            else if(cmpstring(mode,"PRES")|| cmpstring(mode,"RESI"))
            {
                if(cmpstring(token[0],"ATOM"))
                {
                    if( token.size() < 3 ){
                        printf("## SYNTAX ERROR: at least 3 identifiers expected in RESI section, after ATOM: \n");
                        printf("## %s",rawlinestring.c_str());
                        continue;
                    }

                    int Type = findAtomType(token[2]);
                    if(Type<0)
                    {
                        printf("##SYNTAX ERROR: Unknown AtomType: '%s' \n",token[2].c_str());
                        printf("## %s",rawlinestring.c_str());
                        return -1;
                    }

                    // The atom paraemter inherits all of it's Type's parameters
                    AtomParameter newatom(AtomType[Type]);

                    newatom.rawname = token[1];
                    newatom.pdbname = token[1];
                    newatom.type_name = token[2];
                    newatom.comment = comment;
                    newatom.FFType = Type;

                    newatom.valid = 1;

                    double charge;

                    charge = 0.0;
                    if(str2double(token[3],charge)!=0)
                    {
                        printf("## SYNTAX ERROR: Expected a charge (numerical value) as 3rd parameter: \n");
                        printf("## %s",rawlinestring.c_str());
                    }
                    newatom.charge = charge;

                    newmoldef.atom.push_back(newatom);
                    continue;
                }
                else if(cmpstring(token[0],"BOND"))
                {
                    for(int dset=3;dset<=13;dset+=2)
                    {
                        if( token.size() >= dset )
                        {
                            std::string bondnamei = token[dset-2];
                            std::string bondnamej = token[dset-1];

                            int bondi = newmoldef.findAtomRaw(bondnamei);
                            int bondj = newmoldef.findAtomRaw(bondnamej);

                            if(bondnamei.c_str()[0]=='+')
                            {
                                bondi = newmoldef.findAtomRaw(bondnamei.substr(1));
                                if((bondnamej.c_str()[0]!='+') && (bondnamej.c_str()[0]!='-'))
                                    bondnamej = "-" + bondnamej;
                            }
                            if(bondnamei.c_str()[0]=='-')
                            {
                                bondi = newmoldef.findAtomRaw(bondnamei.substr(1));
                                if((bondnamej.c_str()[0]!='+') && (bondnamej.c_str()[0]!='-'))
                                    bondnamej = "+" + bondnamej;
                            }
                            if(bondnamej.c_str()[0]=='+')
                            {
                                bondj = newmoldef.findAtomRaw(bondnamej.substr(1));
                                if((bondnamei.c_str()[0]!='+') && (bondnamei.c_str()[0]!='-'))
                                    bondnamei = "-" + bondnamei;
                            }
                            if(bondnamej.c_str()[0]=='-')
                            {
                                bondj = newmoldef.findAtomRaw(bondnamej.substr(1));
                                if((bondnamei.c_str()[0]!='+') && (bondnamei.c_str()[0]!='-'))
                                    bondnamei = "+" + bondnamei;
                            }
                            if(bondi<0)
                            {
                                printf("##SYNTAX ERROR: Unknown Atom Name '%s' \n",bondnamei.c_str());
                                printf("## %s",rawlinestring.c_str());
                                continue;
                            }
                            if(bondj<0)
                            {
                                printf("##SYNTAX ERROR: Unknown Atom Name '%s' \n",bondnamej.c_str());
                                printf("## %s",rawlinestring.c_str());
                                continue;
                            }
                            if(bondi>=0)
                            {
                                if(!cmpstring(newmoldef.atom[bondi].restraint,""))
                                    newmoldef.atom[bondi].restraint += ",";
                                newmoldef.atom[bondi].restraint += bondnamej;
                            }
                            if(bondj>=0)
                            {
                                if(!cmpstring(newmoldef.atom[bondj].restraint,""))
                                    newmoldef.atom[bondj].restraint += ",";
                                newmoldef.atom[bondj].restraint += bondnamei;
                            }
                        }
                    }
                    continue;
                }
                else if(cmpstring(command,"IMPR") || cmpstring(command,"IMPROPER")) 
                {
                    if( token.size() < 5 )
                    {
                        printf("## SYNTAX ERROR: Expected 4 atom Type identifiers after IMPR \n");
                        printf("## %s",rawlinestring.c_str());
                        continue;
                    }
                    if( token.size() >= 5 )
                    {
                        DihedralDefinition newimpropertype;
                        std::string improperline;
                        improperline = "";
                        improperline += "IMPROPER ";
                        improperline += token[1] + " ";
                        improperline += token[2] + " ";
                        improperline += token[3] + " ";
                        improperline += token[4] + " ";
                        if(newimpropertype.readDefinitionLine(improperline, errlog) == 0) 
                        {
                            newmoldef.improper.push_back(newimpropertype);
                        }
                        else
                        {
                            printf("## ERROR in improper line \n");
                        }
                    }
                    if( token.size() >= 9 )
                    {
                        DihedralDefinition newimpropertype;
                        std::string improperline;
                        improperline = "";
                        improperline += "IMPROPER ";
                        improperline += token[5] + " ";
                        improperline += token[6] + " ";
                        improperline += token[7] + " ";
                        improperline += token[8] + " ";
                        if(newimpropertype.readDefinitionLine(improperline, errlog) == 0) 
                        {
                            newmoldef.improper.push_back(newimpropertype);
                        }
                        else
                        {
                            printf("## ERROR in improper line \n");
                        }
                    }
                }
                else if(cmpstring(token[0],"GROUP"))
                {                    
                    printf("## Ignoring 'GROUP' %s ", linestring.c_str());
                    continue;
                }
                else if(cmpstring(token[0],"DOUBLE"))
                {
                    printf("## Ignoring 'DOUBLE' %s ", linestring.c_str());
                    continue;
                }
                else if(cmpstring(token[0],"TRIPLE"))
                {
                    printf("## Ignoring 'TRIPLE' %s ", linestring.c_str());
                    continue;
                }
                else if(cmpstring(token[0],"ACCE") || cmpstring(token[0],"ACCEPTOR"))
                {
                    printf("## Ignoring 'ACCE' %s ", linestring.c_str());
                    continue;
                }
                else if(cmpstring(token[0],"DONO") || cmpstring(token[0],"DONOR"))
                {
                    printf("## Ignoring 'DONO' %s ", linestring.c_str());
                    continue;
                } 
                else if(cmpstring(token[0],"IC"))
                {
                    printf("## Ignoring 'IC' %s ", linestring.c_str());
                    continue;
                } 
                else
                {
                    printf("## SYNTAX ERROR: Dont understand identifier: '%s' \n", token[0].c_str());
                    printf("## %s",rawlinestring.c_str());
                }
            }
            else
            {
                printf("## CHARMM Importer doesnt understand: \'\n");
                printf("## %s",rawlinestring.c_str());

            }
        }
        if( comment.length() > 0 )
            printf("# %s \n",comment.c_str());

        /*
        if(cmpstring(command,"ALIAS")) { // ALIAS
        if(token.size() < 3) {
        printf("SYNTAX ERROR: Two identifiers expected after 'ALIAS' \n");
        errorstatus=-1;
        continue;
        }
        ResidueAliasDefinition newalias; // save this alias
        newalias.set(token[1],token[2]);
        AliasDef.push_back(newalias);

        } else if(cmpstring(command,"NOAUTOBONDS")) { // AUTOANGLES - instruct to autogenerate angles
        autobonds = false;
        } else if(cmpstring(command,"NOAUTOANGLES")) { // AUTOANGLES - instruct to autogenerate angles
        autoangles = false;
        } else if(cmpstring(command,"NOAUTOTORSIONS")) { // AUTOTORSIONS - instruct to autogenerate torsions
        autotorsions = false;
        } else if(cmpstring(command,"AUTOBONDS")) { // AUTOANGLES - instruct to autogenerate angles
        autobonds = true;
        } else if(cmpstring(command,"AUTOANGLES")) { // AUTOANGLES - instruct to autogenerate angles
        autoangles = true;
        } else if(cmpstring(command,"AUTOTORSIONS")) { // AUTOTORSIONS - instruct to autogenerate torsions
        autotorsions = true;
        } else if(cmpstring(command,"VDW14SCALING")) {
        if(token.size() < 2) {
        printf("SYNTAX ERROR: Identifier expected after '%s' \n", command.c_str());
        errorstatus=-1;
        continue;
        }

        if(str2double(token[1],value)!=0){
        printf("SYNTAX ERROR: Float expected after '%s' \n", command.c_str());
        errorstatus=-1;
        continue;
        }
        Vdw14Scaling = value;
        } else if(cmpstring(command,"ELEC14SCALING")) {
        if(token.size() < 2) {
        printf("SYNTAX ERROR: Identifier expected after '%s' \n", command.c_str());
        errorstatus=-1;
        continue;
        }

        if(str2double(token[1],value)!=0){
        printf("SYNTAX ERROR: Float expected after '%s' \n", command.c_str());
        errorstatus=-1;
        continue;
        }
        Elec14Scaling = value;
        } else if(cmpstring(command,"SECTION")) { // ADDITIONAL SECTION syntax reading/checking
        if(token.size() < 2) {
        printf("SYNTAX ERROR: Identifier expected after '%s' \n", command.c_str());
        errorstatus=-1;
        continue;
        }

        Section newsection(token[1],filename,line);

        newsection.readSectionBlock(file,line);
        section.push_back(newsection);

        } else if(cmpstring(command,"TYPE")) { // TYPE syntax reading/checking
        AtomTypeParameter newatomtype;

        if(newatomtype.AtomTypeReadLine(linestring.c_str(), line)!=0){
        printf("ERROR: Error reading Type entry (line %d) \n",line);
        errorstatus=-1;
        }else{
        AtomType.push_back(newatomtype);
        }

        } else if(cmpstring(command,"BOND")) { // BOND syntax reading/checking
        BondTypeParameter newbondtype;
        if(newbondtype.readDefinitionLine(linestring.c_str(), line, this) == 0) {
        BondType.push_back(newbondtype);
        }else{
        printf("ERROR: Error reading bond entry (line %d) \n",line);
        errorstatus=-1;
        }

        } else if(cmpstring(command,"ANGLE")) { // ANGLE syntax reading/checking
        AngleTypeParameter newangletype;
        if(newangletype.readDefinitionLine(linestring.c_str(), line, this) == 0) {
        AngleType.push_back(newangletype);
        }else{
        printf("ERROR: Error reading angle entry (line %d) \n",line);
        errorstatus=-1;
        }

        } else if((cmpstring(command,"TORSION")) || (cmpstring(command,"IMPROPER"))) { // TORSION/IMPROPER syntax reading/checking
        TorsionTypeParameter newtorsiontype;
        int treturn = newtorsiontype.readDefinitionLine(TorsionType, TorsionType.size(), linestring, line, this);
        if(treturn==0){
        TorsionType.push_back(newtorsiontype);
        }else if (treturn<0){
        printf("ERROR: Error reading torsion entry (line %d) \n",line);
        errorstatus=-1;
        }
        } else if(cmpstring(command,"MOLECULE")) { // MOLECULE syntax reading/checking
        if(token.size() < 2) {
        printf("SYNTAX ERROR: Identifier expected after '%s' \n", command.c_str());
        errorstatus=-1;
        continue;
        }

        MoleculeDefinition newmoleculetype;
        if(newmoleculetype.readMoleculeBlock(file, line, this) != 0) { // read atom information (ENDMOLECULE terminates)
        errorstatus = -1;
        } else {

        // either only name, or a letter followed by a name is given

        if(token.size() == 2){
        newmoleculetype.name = token[1]; // add name
        } else {
        newmoleculetype.name = token[2]; // add name
        newmoleculetype.letter = token[1].c_str()[0]; // add letter
        }

        // create 3 letter name pointer
        int alen;
        alen = (int) strlen(&newmoleculetype.name[0]) - 3;
        if(alen < 0)
        alen = 0;
        newmoleculetype.l3_name = &newmoleculetype.name[alen];

        molecule.push_back(newmoleculetype);

        }

        } else {
        printf("SYNTAX ERROR: Unknown identifier '%s'(line %d)\n", command.c_str(), line);
        errorstatus=-1;
        }

        */
    }

    fclose(file);
    return errorstatus;
}


int FFParamSet::readCHARMMrtf(const std::string &filename)
{
    return 0;
}


