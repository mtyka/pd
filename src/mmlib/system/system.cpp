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

// Precompiled Header
#include "global.h"

// MMLib Headers
#include "maths/maths.h"
#include "forcefields/ffparam.h"
#include "forcefields/forcefield.h"
#include "workspace/space.h"
#include "fileio/pdb.h"

// Own Header
#include "system/system.h"

// Namespace Declarations
using namespace Maths;

System::System(): ptr_ffps(&nullffps)
{
}

System::System(const FFParamSet &_ptr_ffps):
	ptr_ffps(&_ptr_ffps)
{
};

System::~System()
{
}

Molecule &System::getMolecule(size_t index)
{
	return m_Molecule[index];
}

const Molecule &System::getMolecule(size_t index) const
{
	return m_Molecule[index];
}

size_t System::nMolecules() const
{
	return m_Molecule.size(); 
}
size_t System::nAtoms() const
{
	size_t _natoms = 0;
	for(size_t imol=0;imol<nMolecules();imol++) _natoms += getMolecule(imol).nAtoms();
  return _natoms; 
}

/// \details before adding a new MoleculeBase:
/// First check if we have any molecules loaded - if not then the 
/// molecule dictates what ffps this system will connect to, i.e. 
/// the private ptr_ffps is set to that of the molecule
/// If Molecules have already been loaded then compare the ffps of
/// the molecule to the one of the system - they must match 
/// (i.e. be identical objects with identical addresses !!)
/// Otherwise throw an exceptions
/// \author M.Tyka
void System::add(const Molecule &newmol)
{
	if(m_Molecule.size() == 0){
		ptr_ffps = &newmol.ffps(); // set System::ptr_ffps to that of molecule
	}else{
		if( &newmol.ffps() != ptr_ffps ){
			throw(ArgumentException("\
Cannot add molecule to system because it was created with a different FFParamSet\n\
Molecules belonging to different parameter sets cannot be mixed in one simulation\n\
system for understandable reasons. \n"));
		}
	}
	m_Molecule.push_back(newmol);
}

void System::add(const System &newSys)
{
	if(m_Molecule.size() == 0){
		ptr_ffps = &newSys.ffps(); // set System::ptr_ffps to that of molecule
	}else{
		if( &newSys.ffps() != ptr_ffps ){
			throw(ArgumentException("\
Cannot add system to system because it was created with a different FFParamSet\n\
Molecules belonging to different parameter sets cannot be mixed in one simulation\n\
system for understandable reasons. \n"));
		}
	}
	for( size_t i = 0; i < newSys.nMolecules(); i++ )
	{
		m_Molecule.push_back(newSys.getMolecule(i));
	}
}


/// Adds a Particle to this System as an individual MoleculeBase (it creates a copy)
void System::add(const Particle &newatom, const FFParamSet &ffps)
{
	Molecule newmol(ffps);
	newmol.addParticle(newatom);
	m_Molecule.push_back(newmol);
}



void System::remove(size_t index){
	if(index < m_Molecule.size() )
	{
		m_Molecule.erase( m_Molecule.begin()+ index );
	}
	else
	{
		throw(ArgumentException("Molecule Index out of range\n"));
	}
}




/// Adds lots of copies of copies of molecule to this system filling
/// up the space layed out by boundary
void System::solvate_N(
	const Molecule &newmol,
	const ClosedSpace &boundary, 
	unsigned N,
	int      enforceN
)
{

	if(enforceN > ((float)N*1.5) ){
		printf("Error: enforceN must not be more than 150 percent of N when creating a solvent box\n");
		printf("You asked for certain density number (N = %d) of solvent molecules but asked for *exactly* \n(enforceN = %d ) of solvent molecules to be added.\n",N,enforceN);
		throw(ProcedureException("Error occured during solvent box creation (see error message above)"));
	}
	std::vector<Maths::dvector> points;
	boundary.getEvenPointDistribution(points, N);

	size_t i;
	Molecule tempmol(newmol);
	tempmol.zeroCentreOfGeometry();
	printf("Creating solvent box using a simulation box: %s \n", newmol.name.c_str());
	boundary.info();

	double cutoff=3.0;
	size_t imol_original = nMolecules();

	unsigned countSolventMols = 0;
	for(i=0;i<points.size();i++){
		// check the point is not obscured

		bool reject = false;
		for(size_t imol=0;imol<imol_original;imol++){
			for(int iat = 0; iat < getMolecule(imol).nAtoms(); iat++) {
				//if(dist(getMolecule(imol).atomxyz(iat), points[i]) < cutoff){
				if(dist(getMolecule(imol).atomxyz(iat), points[i]) < /*cutoff*/ max( cutoff, getMolecule(imol).getAtom(iat).radius ) ){
					reject = true;
					break;
				}
			}
			if(reject) break;
		}
		if(reject) continue;
		Molecule tempmol2(tempmol);
		tempmol2.move(points[i]);
		matrix3x3 rmat;
		rmat.setToRandomRot();
		tempmol2.rotate(rmat);
		countSolventMols ++;
		add(tempmol2);
		
		// If user specified the exact number of solvent molecules and we're over then stop adding.
		if ( enforceN >= 0 )
		if( countSolventMols >= enforceN) break;
	}

	// If user specified the exact number of solvent molecules and we havnt added enough then add more :D
	if ( enforceN >= 0 )
  while( countSolventMols < enforceN){
    // go through list randomly and insert more solvent molecules at points intermediate to the
    // original positions (to minimise clashes) 
		i = rand()%(points.size()-2) + 1;
		dvector insertvectorb = points[i-1];
		dvector insertvector = points[i];
		dvector insertvectora = points[i+1];
		if(dist(insertvector,insertvectora) < dist(insertvectorb,insertvector)){
			insertvector.add(insertvectora);
		}else{
			insertvector.add(insertvectorb);
		}		
		insertvector.mul(0.5);

		bool reject = false;
		for(size_t imol=0;imol<imol_original;imol++){
			for(int iat = 0; iat < getMolecule(imol).nAtoms(); iat++) {
				
				if(dist(getMolecule(imol).atomxyz(iat), insertvector) < /*cutoff*/ max( cutoff, getMolecule(imol).getAtom(iat).radius ) ){
					reject = true;
					break;
				}
			}
			if(reject) break;
		}
		if(reject) continue;
		Molecule tempmol2(tempmol);
		tempmol2.move(insertvector);
		matrix3x3 rmat;
		rmat.setToRandomRot();
		tempmol2.rotate(rmat);
		countSolventMols ++;
		add(tempmol2);
	}

	printf("  Created solvent box: %d molecules added to system\n",countSolventMols);

  
}

void System::solvate(
	const Molecule &newmol, 
	const ClosedSpace &boundary, 
	double density,               // in g/ml, e.g. water = 1.0
	int    enforceN
)
{
	double volume = boundary.volume() * cube(Physics::PhysicsConst::Angstrom); // in m^3
	double M =  1.0 / (newmol.getTotalMass()); // #/kg
	double dblN = unsigned( density * 1000.0 // kg/m^3
													* M * volume );
	unsigned N = (unsigned) dblN;

	solvate_N(newmol, boundary, N, enforceN);
}

int System::setup() 
{
	if(checkAllMoleculeParams()<0){
		return loadAllMoleculeParams();
	}
	return 0;
}

void System::info() const 
{
	printf("System Info: --- %s -------------------------------\n","System"); //name.c_str());
	unsigned total_atoms = 0;
	unsigned total_res = 0;
	unsigned total_mol = 0;
	double   total_mass = 0;
	double   total_charge = 0;
	printf("             No  .----- Name -----.    Atoms   Res    Mols     Mass    Charge\n"); 
	int identical_count=1;
	char buffer[1024];
	for(size_t i=0;i<m_Molecule.size();i++)
	{
		if(i<(m_Molecule.size() - 1)){
			if( m_Molecule[i] == m_Molecule[i+1] ){
				identical_count += 1;
				continue;
			}
		}
		if(identical_count <= 1){
			sprintf(&buffer[0],"%d",i);
		}else{
			sprintf(&buffer[0],"%d-%d",i-identical_count+1,i);
		}
		printf("%9s ",&buffer[0]);
		printf("%5d%20s %8d %6d %4d %8.1lf  %+8.2lf \n", 
			identical_count,
			m_Molecule[i].name.c_str(),
			m_Molecule[i].nAtoms(),
			m_Molecule[i].nResidues(),
			m_Molecule[i].nMolecules(),
			m_Molecule[i].getTotalMass() / Physics::PhysicsConst::amu,
			m_Molecule[i].getTotalCharge() );

		total_atoms  += m_Molecule[i].nAtoms()*identical_count;
		total_res    += m_Molecule[i].nResidues()*identical_count;
		total_mol    += m_Molecule[i].nMolecules()*identical_count;
		total_mass   += m_Molecule[i].getTotalMass() / Physics::PhysicsConst::amu*identical_count;
		total_charge +=	m_Molecule[i].getTotalCharge()*identical_count;

		identical_count = 1;
	}
	printf("-------------------------------------------------------------------------\n"); 
	printf("           %4d                     %8d %6d %4d %8.1lf  %+8.2lf \n", 
			m_Molecule.size(),
			total_atoms,  
			total_res,    
			total_mol,    
			total_mass,   
			total_charge); 


}

void System::detail() const
{
	printf("System Details: --------- \n");
	for(size_t i=0;i<m_Molecule.size();i++)
	{
		m_Molecule[i].detail();
	}
}

int System::loadAllMoleculeParams()
{
	printf("Loading all system molecule parameters..\n");
	for(size_t imol=0;imol<nMolecules();imol++)
	{
		if(getMolecule(imol).loadParameters(ffps())!=0)
		{
			printf("ERROR: Error loading parameters for molecule %d\n",imol);
			return -1;
		}
	}
	printf("Done loading .. \n");
	return 0;
}

int System::checkAllMoleculeParams()
{
	printf("Checking all system molecule parameters.. %d molecules\n",nMolecules());
	for(size_t imol=0;imol<nMolecules();imol++)
	{
		if(getMolecule(imol).checkParameters()!=0)
		{
			return -1;
		}
	}
	printf("Done checking ..\n");
	return 0;
}

// gets maximum absolute coordinate in system (i.e. largest value of |x|,|y| and |z|
dvector System::getEncompassingVector() const
{
	dvector largest(0,0,0);
	for(size_t imol=0;imol<nMolecules();imol++)
	{
		dvector vec = getMolecule(imol).getEncompassingVector();
		if(fabs(vec.x) > largest.x ) largest.x = fabs(vec.x);
		if(fabs(vec.y) > largest.y ) largest.y = fabs(vec.y);
		if(fabs(vec.z) > largest.z ) largest.z = fabs(vec.z);
	}
	return largest;
}

double System::getTotalMass() const
{
	double totalmass = 0;
	for(size_t imol=0;imol<nMolecules();imol++){
		totalmass += getMolecule(imol).getTotalMass(); 
	}
	return totalmass;
}

double System::getTotalCharge() const
{
	double totalcharge = 0;
	for(size_t imol=0;imol<nMolecules();imol++){
		totalcharge += getMolecule(imol).getTotalCharge();
	}
	return totalcharge;
}

dvector System::getCentreOfMass() const
{
	dvector com(0,0,0);
	dvector p;
	for(size_t imol=0;imol<nMolecules();imol++){
		for(int i = 0; i < getMolecule(imol).nAtoms(); i++) {
			p.setTo(getMolecule(imol).atomxyz(i));
			p.mul(getMolecule(imol).atom[i].mass);
			com.add(p);
		}
	}
	com.div(getTotalMass());
	return com;
}

dvector System::getCentreOfGeometry() const
{
	dvector cog(0,0,0);
	for(size_t imol=0;imol<nMolecules();imol++){
		for(int i = 0; i < getMolecule(imol).nAtoms(); i++){
			cog.add(getMolecule(imol).atomxyz(i));
		}
	}
	
	cog.div((double)nAtoms());
	return cog;
}




void System::calcInertiaTensor(matrix3x3 & I) const
{
	matrix3x3 Ii;
	I.setToNull();
	dvector posSI;
	dvector com = getCentreOfMass();
	for(size_t imol=0;imol<nMolecules();imol++){
		for(int i = 0; i < getMolecule(imol).nAtoms(); i++) {
			posSI.setTo(getMolecule(imol).atomxyz(i));
			posSI.sub(com);
			posSI.mul(Physics::PhysicsConst::Angstrom);
			Ii.setTo( sqr(posSI.y) + sqr(posSI.z), -posSI.x * posSI.y,						-posSI.x * posSI.z,
							 -posSI.y * posSI.x					 ,  sqr(posSI.x) + sqr(posSI.z),	-posSI.y * posSI.z,
							 -posSI.z * posSI.x					 , -posSI.z * posSI.y,						 sqr(posSI.x) + sqr(posSI.y));
			Ii.mul(getMolecule(imol).atom[i].mass);
			I.add(Ii);
		}
	}
}

	
void System::alignAlongPrincipalAxes() 
{
	matrix3x3 I;
	calcInertiaTensor(I);
	double lambda1,lambda2,lambda3;
	dvector p1,p2,p3;
	I.diagonaliseSymetric(lambda1,lambda2,lambda3,p1,p2,p3);
	matrix3x3 rot;
	dvector i(1,0,0);
  dvector j(0,1,0);
	// rotate the system so that the first tow principal 
	// components point along the cartesian axe
  // by definition the last will also point inthe right direction
	superimpose(i,j,p1,p2,rot);   
	rotate(rot);
}

void System::printInertiaInfo() const
{
	matrix3x3 I;
  calcInertiaTensor(I);
	printf("Ixx =        %.6e   kg.m^2 \n", I.r[0][0] ); 
	printf("Iyy =        %.6e   kg.m^2 \n", I.r[1][1] ); 
	printf("Izz =        %.6e   kg.m^2 \n", I.r[2][2] ); 
	printf("Ixy = Iyx =  %.6e   kg.m^2 \n", I.r[0][1] ); 
	printf("Ixz = Izx =  %.6e   kg.m^2 \n", I.r[0][2] ); 
	printf("Izy = Iyz =  %.6e   kg.m^2 \n", I.r[2][1] ); 
  printf("I = IA IB IC =  %.6e  (kg.m^2)^3 \n",
           I.r[0][0]*I.r[1][1]*I.r[2][2] +
           2.0 * I.r[0][1]*I.r[0][2]*I.r[2][1] -
           I.r[0][0] * sqr(I.r[2][1]) -
           I.r[1][1] * sqr(I.r[2][0]) -
           I.r[2][2] * sqr(I.r[0][1]) );
} 


double System::calcRotationalPartition(double temp, unsigned SymNumber) const
{
  matrix3x3 I;
  calcInertiaTensor(I);

  // product of principal moments of inertia
  double IAIBIC =
           I.r[0][0]*I.r[1][1]*I.r[2][2] +
           2.0 * I.r[0][1]*I.r[0][2]*I.r[2][1] -
           I.r[0][0] * sqr(I.r[2][1]) -
           I.r[1][1] * sqr(I.r[2][0]) -
           I.r[2][2] * sqr(I.r[0][1]);
	double qrot;
	qrot = 8.0*sqr(MathConst::PI)/(double(SymNumber));
	qrot *= sqrt(cube(2.0*MathConst::PI * Physics::PhysicsConst::kB * temp/sqr(Physics::PhysicsConst::planck))*IAIBIC);
	return qrot;
}


void System::zeroCentreOfGeometry(){
	int i;
	dvector cog = getCentreOfGeometry();
	for(size_t imol=0;imol<nMolecules();imol++){
		for(i = 0; i < getMolecule(imol).nAtoms(); i++) {
			getMolecule(imol).atomxyz(i).sub(cog);
		}
	}

}

void System::zeroCentreOfMass(){
	int i;
	dvector com = getCentreOfMass();
	for(size_t imol=0;imol<nMolecules();imol++){
		for(i = 0; i < getMolecule(imol).nAtoms(); i++) {
			getMolecule(imol).atomxyz(i).sub(com);
		}
	}
}

void System::rotate(const matrix3x3 &rot)
{
	for(size_t imol=0;imol<nMolecules();imol++){
		for(int i = 0; i < getMolecule(imol).nAtoms(); i++) {
			getMolecule(imol).atomxyz(i).mulmat(rot);
		}
	}
}


void System::save( IO::OutputFile &_output ){
	_output.save( *this );
}


////////////////////////////////////////////////////////////////////////////////
//
// Debug Tools & Short Cuts. 

void System::printPDB(const std::string& _Filename)
{
	IO::PDB_Writer out( _Filename, false );
	out.write(*this);
}

void System::printPDB()
{
	IO::PDB_Writer out( std::cout, false );
	out.write(*this);
}

