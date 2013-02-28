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

#ifndef __SYSTEM_H
#define __SYSTEM_H

// System Includes
#include <stdlib.h>
#include <string>
#include <vector>

// MMLib Includes
#include "forcefields/ffparam.h"
#include "system/fundamentals.h"
#include "fileio/outtra.h"
#include "system/molecule.h" // Provides member

class ClosedSpace;

//-------------------------------------------------
//
/// \brief A System is a container for Molecules and provides functions for manipulating them. 
///
/// \details 
/// A simulation system is set up in a System which is mainly a container that holds objects of the 
/// type Molecule, which in turn possess lists of atom coordinates and properties. The system is thus 
/// stored in a hierarchical manner allowing easy modification and duplication of molecules and other 
/// components. The Molecules are created as separate objects and the topology information (such as 
/// covaleing information and other geometrical parameters) are supplied by the FFParamSet 
/// class. The created molecules are then inserted into the system to create the simulation set up. For 
/// example a protein and a ligand can be loaded, brought into correct relative orientation and finally 
/// solvated with a box of water, as well as a few counter ions. For purpose of simulation this system 
/// is then exported into a WorkSpace in which the atoms are stored in a linear array for faster and 
/// simpler access by the Forcefields and Protocols.
///
/// \author Mike Tyka & Jon Rea 
///
class PD_API System 
{
public:
    System();

    /// Constructor must take a reference to an ffps
    System(const FFParamSet &_ptr_ffps);

    ~System();

    /// returns a non-const reference to a molecules stored
    Molecule &getMolecule(size_t index);

    /// returns a const reference to a molecules stored
    const Molecule &getMolecule(size_t index) const;

    /// returns the number of molecules stored
    size_t nMolecules() const;

    /// returns the total number of atoms 
    size_t nAtoms() const;

    /// Adds a Molecule to this System (it creates a copy)
    void add(const Molecule &newmol);

    /// Adds the contents of another System to this one  
    void add(const System &newSys);

    /// Remove the Molecule with the index 
    void remove(size_t index);

    /// Adds a Particle to this System as an individual MoleculeBase (it creates a copy)
    ///
    /// A FFParamSet must be supplied here because particles dont have their
    /// own reference to a parameter set.
    void add(const Particle &newatom, const FFParamSet &ffps);

    /// Adds lots of copies of copies of molecule to this system filling
    /// up the space layed out by boundary
    void solvate_N(
        const Molecule &newmol, 
        const ClosedSpace &boundary, 
        unsigned N, 
        int enforceN = -1 ///< if enforceN >= 0, Force to insert this number of solvent molecules (useful if you want to create two systems with identical numbers of water atoms
        );

    void solvate(
        const Molecule &newmol, 
        const ClosedSpace &boundary, 
        double density, ///< required density of the solvent in g/ml, e.g. water = 1.0
        int enforceN = -1 ///< if enforceN >= 0, Force to insert this number of solvent molecules (useful if you want to create two systems with identical numbers of water atoms
        );

    /// this can be used to access ffps through the system
    const FFParamSet &ffps() const 
    {
        return *ptr_ffps;
    }

    /// gets maximum absolute coordinate in system (i.e. largest value of |x|, |y| and |z|
    Maths::dvector getEncompassingVector() const;

    /// Returns the total mass of the system
    double getTotalMass() const; 

    /// Returns the total charge of the system
    double getTotalCharge() const; 

    /// Returns the center of mass of the system
    Maths::dvector getCentreOfMass() const; 

    /// Returns the center of geometry of the system
    Maths::dvector getCentreOfGeometry() const; 

    /// calculates the inertia tensor of the entire system
    void calcInertiaTensor(Maths::matrix3x3 & I) const; 

    /// aligns the system such that the principal axes of rotation point along the coordinate axes.
    /// In other words, by applying this function, the smallest solvent box is needed to accomodate the molecule.
    /// Note that this has to be called before the solvent is added. 
    void alignAlongPrincipalAxes();

    /// displays the inertia tensor 
    void printInertiaInfo() const;

    /// calculates the rotational partition function of the system 
    double calcRotationalPartition(double temp, unsigned SymNumber) const;

    /// rotates the entire system by the rotation matrix
    void rotate(const Maths::matrix3x3 &rot);

    /// center the system such that the centre of geometry is a s 0,0,0
    void zeroCentreOfGeometry();

    /// center the system such that the centre of mass is a s 0,0,0
    void zeroCentreOfMass();	

    /// Displays a short summary of the system 
    void info() const;

    /// Displays detailed information about the system
    void detail() const;

    // ----------------
    // IO functionality

    /// This saves the current state of the system in a file format of the users choice
    ///
    ///  For example  mysystem.save( OutputFile_PDB("mypdbfile") ) 
    ///  saves a PDB file.
    void save( IO::OutputFile &_output );

    /// These are dirty shortcuts to dump the state of the System in a PDB file!
    void printPDB(const std::string& _Filename);
    void printPDB();

    friend class WorkspaceCreatorBase;

private:
    int setup();
    int checkAllMoleculeParams();
    int loadAllMoleculeParams();

    std::vector<Molecule> m_Molecule;
    const FFParamSet *ptr_ffps;
};

#endif

