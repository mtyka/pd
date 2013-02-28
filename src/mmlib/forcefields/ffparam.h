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

// ffparam.h
// * ------------------------------------------------------------------------
// Atomic Parameter Loading
// ------------------------------------------------------------------------ */

#ifndef __FFPARAM_H
#define __FFPARAM_H

#include <string>
#include <vector>

#include "system/fundamentals.h"
#include "library/mapper.h"  // required as it provides a base class
#include "library/nameconv.h"  // required as it provides a member variable

#ifndef SWIG

class PD_API FFParamSet; // the main class PD_API structure, defined/declared later
class PD_API MoleculeDefinition;
class PD_API AtomParameter;

/// These are empty pre-allocated objects of the above classes 
/// Used in situations where an "empty" reference is needed.
extern const FFParamSet nullffps;
extern const MoleculeDefinition nullmolecule;
extern const AtomParameter nullatomparameter;


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
class PD_API CovalentLink
{
public:
    CovalentLink()
    {
        i = 0;
        ani = "";
        roi = 0;
    }
    int i; // internal molecule atom index of bond partner
    std::string ani; // name of bonded atom
    int roi; // residue offset (is != 0 if linking to next/previous residue)
};


//-------------------------------------------------
//
/// \brief  used to save atom Type definitions
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
class PD_API AtomTypeParameter: public Particle
{ 
public:
    AtomTypeParameter()
    {
        init();
    }

    /// Name of Atom Type (could the inherited type_name be used here instead?)
    std::string name; 
    bool used;

    int AtomTypeReadLine(const std::string linestr, ParseErrorLogger &errlog);

private:
    void init()
    {
        name = "";
        used = false;
    }
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
class PD_API AtomParameter: public AtomTypeParameter
{
public:
    AtomParameter(): AtomTypeParameter() 
    {
        init();
    }

    AtomParameter( const AtomTypeParameter &copy)
        : AtomTypeParameter(copy) 
    {
        init();
    }

    friend class PD_API FFParamSet;
    friend class PD_API MoleculeDefinition;

    char valid; // 0 for invalid (to be removed/ignored) !=0 for valid	

    std::vector<CovalentLink> r_covalent; // bond definitions

private:
    void init()
    {
        comment = ""; // comment
        valid = 1;
    }

    std::string restraint;
    std::string comment;
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
class PD_API BondTypeParameter
{
public:

    /// Internal type number of second atom type
    int i; 

    /// Internal type number of second atom type
    int j; 
    char used;

    /// Equilibrium bond length
    double length;    

    /// forceconstant for harmonic bond restraints
    double forceconstant;  

    /// order of bond - 1 = single, 2 = double, 3 = triple
    double bondorder;     

    /// Reads a line frm the forcefield definition file that contains a bond type entry:
    /// For example "BOND  CA  N   1.35  300  1"
    int readDefinitionLine(const std::string linestring, FFParamSet &ffps, ParseErrorLogger &errlog);

    /// compares if the bond Type matches t{i,j}
    int cmpBond(int ti, int tj) const; 
};


//-------------------------------------------------
//
/// \brief Stores information about a type of angle 
///
/// \details It holds three integers holding the indices of the atom types concerned (stored in an FFParamSet) 
///   as well as the equilibrium angle and the forceconstant    
/// \author Mike Tyka & Jon Rea 
///
class PD_API AngleTypeParameter
{
public:
    char used;

    /// Type identifier (in AtomType[])
    int i, a, j; 						
    double angle;
    double forceconstant;

    int readDefinitionLine(const std::string linestring, FFParamSet &ffps, ParseErrorLogger &errlog);

    /// compares if the angle yype matches t{i,a,j}
    int cmpAngle(int ti, int ta, int tj) const;
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
class PD_API TorsionTypeParameter
{
public:
    TorsionTypeParameter()
    {
        zero();
    }

    char used;
    int i, a, b, j;

    /// Type identifier (in AtomType[]) 0 for torsion, 1 for improper
    int Type;

    /// how many fourier terms
    int terms; 

    /// the individual parameters for each term
    double Vn[4]; 
    double n[4];
    double gamma[4];

    int addTerm(double newVn, double newn, double newgamma);
    int cmpTorsion(int ti, int ta, int tb, int tj, int Type) const;
    void zero();
    void info();

    int readDefinitionLine(
        std::vector<TorsionTypeParameter> &list,
        int list_length,
        const std::string &linestring,
        FFParamSet &ffps,
        ParseErrorLogger &errlog
        );
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
class PD_API DihedralDefinition
{
public:

    DihedralDefinition()
    {
        roi = 0;
        roa = 0;
        rob = 0;
        roj = 0;
        Type = 0;
    }

    std::string ani; // Atom Names
    std::string ana;
    std::string anb;
    std::string anj;

    int roi;
    int roa;
    int rob;
    int roj; // residue offsets (allows inter residue definitions

    int Type;

    int readDefinitionLine(const std::string &linestring, ParseErrorLogger &errlog);
    int printDefinitionLine(FILE * file);
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
class PD_API MoleculeDefinition
{
public:
    friend class PD_API FFParamSet;

    MoleculeDefinition() 
        : letter('-'),
        backpos(),
        frwdpos(),
        m_AtomAppendMode(true),		
        l3_name("---"),
        name("unknown")
    {
    }

    MoleculeDefinition(const MoleculeDefinition &newmol)
    {
        (*this) = newmol;
    }

    MoleculeDefinition &operator=(const MoleculeDefinition &newmol)
    {		
        m_AtomAppendMode = newmol.m_AtomAppendMode;

        letter = newmol.letter;
        name = newmol.name;
        l3_name = newmol.l3_name;

        backName = newmol.backName;
        frwdName = newmol.frwdName;
        backpos = newmol.backpos;
        frwdpos = newmol.frwdpos;

        atom = newmol.atom;
        improper = newmol.improper;
        torsion = newmol.torsion;

        // change the parent links back to "this"
        for(size_t i=0;i<atom.size();i++) 
            atom[i].parent = this;

        return (*this);
    }

    int findAtomRaw(const std::string &atom_name) const;
    int findAtomPDB(const std::string &atom_name) const;
    int findImproperDef(const DihedralDefinition &query) const;
    int checkRestraintList( ParseErrorLogger &errlog );
    int decodeRestraintList(
        const std::string &data,
        std::vector<CovalentLink> &bondlist,
        ParseErrorLogger &errlog 
        );

    int decodeRestraintList(  ParseErrorLogger &errlog );

    int readMoleculeBlock(FILE * file, ParseErrorLogger &errlog, FFParamSet &ffps);
    int writeMoleculeBlock(FILE *file, FFParamSet &ffps);

    // inspectors
    const char *c_l3_name() const
    {
        return l3_name.c_str();
    }

    const std::string &s_l3_name() const
    {
        return l3_name;
    }

    const char *c_name() const
    {
        return name.c_str();
    }

    const std::string &s_name() const
    {
        return name;
    }

public:
    char letter;

    /// False if no BackLink, true if links to prev. residue
    bool hasBackLink() const;
    /// False if no FrwdLink, true if links to next residue
    bool hasFrwdLink() const;

    std::string backName;
    std::string frwdName;
    Maths::dvector backpos;
    Maths::dvector frwdpos;

    std::vector<AtomParameter> atom; ///< stores Atom information
    std::vector<DihedralDefinition> improper; ///< stores explicit IMPROPER information
    std::vector<DihedralDefinition> torsion; ///< stores explicit TORSION information

private:
    std::string name;
    std::string l3_name;
    std::string comment;

    /// Controlls wether new atoms definitions are added at the start or end of their list. This is because atom order is important for rotations.
    bool m_AtomAppendMode; 
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
class PD_API Section
{
public:
    Section(){}
    Section(const std::string &_sectionname,
        const std::string &_filename,
        int _lineoffset)
    {
        sectionname = _sectionname;
        filename = _filename;
        lineoffset = _lineoffset;
    }
    ~Section(){}

    int readSectionBlock(FILE * file, ParseErrorLogger &errlog);

    std::vector< std::string > sectionline;
    int lineoffset;
    std::string filename;
    std::string sectionname;
};

#endif


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
class PD_API FFParamSet : 
    public Object, 
    public Library::AliasMapper, 
    public Library::ClassMapper 
{
public:
    FFParamSet() 
        : Object(), 
        Library::AliasMapper(), 
        Library::ClassMapper() 
    {
        settodefault();
    }

    FFParamSet(const std::string  &_filename) 
        : Object(), 
        Library::AliasMapper(), 
        Library::ClassMapper() 
    {
        settodefault();
        readLib(_filename);
    }

    virtual FFParamSet* clone() const { return new FFParamSet(*this); }
    virtual ~FFParamSet(){}
    void settodefault();

public:

#ifndef SWIG
    std::string fffilename;
    std::string ffidentifier;

    std::vector<MoleculeDefinition> molecule;
    std::vector<AtomTypeParameter> AtomType;
    std::vector<BondTypeParameter> BondType;
    std::vector<AngleTypeParameter> AngleType;
    std::vector<TorsionTypeParameter> TorsionType;

    std::vector<Section> section;
    std::vector<Library::ResidueAliasDefinition> AliasDef;
    std::vector<std::string> m_CustomProperty;
    std::vector<std::string> m_ReservedKeywords;

    // Various flags & parameters
    bool identgiven;
    bool autobonds;
    bool autoangles;
    bool autotorsions;

    double Vdw14Scaling;
    double Elec14Scaling;
#endif

    int readLib(const std::string &filename, int level = 0);
    int writeLib(const std::string &filename);

    int readCHARMMprm(const std::string &filename);
    int readCHARMMrtf(const std::string &filename);

    int unifyForcefield();

#ifndef SWIG
    void checkUsedTypes();

    int findAtomType(const std::string &name) const;
    int findBondType(int ti, int tj) const;
    int findAngleType(int ti, int ta, int tj) const;
    int findTorsionType(int ti, int ta, int tb, int tj) const;
    int findImproperType(int ti, int ta, int tb, int tj) const;
    int findSection(const std::string &name) const;
    int findCustomProperty(const std::string &name) const;
    int findReservedKeyword(const std::string &name) const;

    int findAliasName(const std::string &_queryAlias,std::string &_returnAliasName) const;
    int findMoleculeType_withoutAlias(const std::string &queryname) const;
    int findMoleculeType(const std::string &MolName) const;

    int getSection(const std::string &queryname, Section &result) const;

#endif

    // AliasMapper implementation

    /// Look up the long name for a given single letter name
    virtual bool lookupLongName( char _SourceID, std::string& _DestID ) const; 

    /// Look up the single letter for a long residue name
    virtual bool lookupShortName( const std::string& _SourceID, char& _DestID ) const; 

    /// Look to see if the source name is an alias, and return the true name
    virtual bool lookupAlias( const std::string& _SourceID, std::string& _DestID ) const; 

private:

    void printAtomParameters(AtomParameter * apar);

};

#endif

