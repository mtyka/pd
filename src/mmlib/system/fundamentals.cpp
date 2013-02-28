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

#include <iomanip>

#include "system/fundamentals.h"
#include "system/molecule.h"
#include "forcefields/ffparam.h"
#include "library/valuestore.h"

const unsigned int Particle::flag_Valid = 1;            ///< All good or something wrong ?
const unsigned int Particle::flag_Used = 2;             ///< Used or not ?
const unsigned int Particle::flag_Dummy = 4;            ///< Real atom or dummy particle ?
const unsigned int Particle::flag_Static = 8;           ///< should it be moved ?
const unsigned int Particle::flag_KnownStructure = 16;  ///< Is this atom position of known structural data ?
const unsigned int Particle::flag_RebuildRequired = 32; ///< Does this atom need rebuilding ?

// Chemical/Structural Nature
const unsigned int Particle::flag_Hydrogen = 256;       ///< Heavy atom or not ?
const unsigned int Particle::flag_Backbone = 512;       ///< backbone or sidechain ?
const unsigned int Particle::flag_C_alpha  = 1024;      ///< Alpha Carbon ?
const unsigned int Particle::flag_CoreBackbone = 2048;  ///< One of the core heavy backbone atoms?

const unsigned int Particle::flag_Moved = 4096;         ///< Moved since the last reset ? 

Particle::Particle()
{
	init();
}

Particle::Particle(const std::string &_rawname,
				   const std::string &_parentname)
{
	init();
	rawname = _rawname;
	parentname = _parentname;
}

const Particle& Particle::operator= ( const Particle& _newpart )
{
	if(_newpart._posPointer == &_newpart._pos){
		_posPointer   =	&_pos;						
	}else{
		_posPointer = _newpart._posPointer;
	}
	_pos         	=	_newpart._pos;         	
	_posRef				=	_newpart._posRef;      
	_posGeom			=	_newpart._posGeom;     
	parent				=	_newpart.parent; 
	i             =	_newpart.i;                      
	ir            =	_newpart.ir;                      
	imol          =	_newpart.imol;                
	igroup        =	_newpart.igroup;                  
	cov12atom			=	_newpart.cov12atom;
	cov13atom			=	_newpart.cov13atom;
	cov14atom			=	_newpart.cov14atom;
	rawname				=	_newpart.rawname;        
	pdbname				=	_newpart.pdbname;        
	type_name			=	_newpart.type_name;      
	parentname		=	_newpart.parentname;     
	parentl3name	=	_newpart.parentl3name;   
	parentletter  =	_newpart.parentletter;           
	Z							=	_newpart.Z;          
	mass					=	_newpart.mass;    
	radius				=	_newpart.radius;  
	epsilon				=	_newpart.epsilon; 
	charge				=	_newpart.charge;  
	epot					=	_newpart.epot;    
	FFType				=	_newpart.FFType;     
	m_CustomProp	=	_newpart.m_CustomProp;
	flags					=	_newpart.flags;

	return (*this);
}

void Particle::CopyParamsFrom(const Particle& _CopyParticle)
{
	Z = _CopyParticle.Z;
	mass = _CopyParticle.mass;
	radius = _CopyParticle.radius;
	epsilon = _CopyParticle.epsilon;
	charge = _CopyParticle.charge;
	FFType = _CopyParticle.FFType;
	parent = _CopyParticle.parent;
}

void Particle::init()
{
	i = 0;
	ir = 0;
	imol = 0;
	igroup = 0;

	flags = 0; // initialise all flags to 0/'off'

	resetPosPointer();  // by default _posPointer will point to the internal position store.

	// set coordinates to 0,0,0
	pos().setTo(0,0,0);
	posRef().setTo(0,0,0);
	posGeom().setTo(0,0,0);
	rawname = "?";
	pdbname = "?";
	type_name = "?";
	parentname = "?";
	parentl3name = "?";
	parentletter = '?';

	Z = 0;
	mass = 0;
	radius = 0;
	epsilon = 0;
	charge = 0;

	epot = 0;
	parent=&nullmolecule;

	// purge
	FFType = 0;
}



void Particle::addProperty(const std::string &_name,const std::string &_data)
{
	Library::ValueStore *vs = Library::ValueStore::getSingleton();

	// check if the property already exists for this atom
	// if this atom is alreaady linked to a property of this name
	// remove the link (but not the property itself since other atoms/things
	// may rely on it
	// Then add a NEW property with that name and link it to this atom
	size_t i_id;
restart:
	for(i_id = 0; i_id < m_CustomProp.size(); i_id ++ )
	{
		if(vs->hasName(m_CustomProp[i_id],_name)){  
			// we found a property under that name
			// get rid of the link
			std::vector< size_t >::iterator iter = m_CustomProp.begin();
			iter += i_id;
			m_CustomProp.erase(iter);
			goto restart; // try again to make sure no more properties under that name remain
		}
	}
	size_t id = vs->addNewData(_name, _data); 
	m_CustomProp.push_back( id );
}

/// Custom Property - obtain the data string
std::string Particle::getCustomProperty_String(const std::string &_name) const {
	Library::ValueStore *vs = Library::ValueStore::getSingleton();

	size_t i_id;
	for(i_id = 0; i_id < m_CustomProp.size(); i_id ++ )
	{
		if(vs->hasName(m_CustomProp[i_id],_name)){
			std::string data = vs->getRawData(m_CustomProp[i_id]);
			return data;
		}
	}

	throw(ProcedureException( std::string("Atom ") + pdbname + 
		std::string(" does not have a custom property ") + 
		_name + 
		std::string(". Check your script and the input forcefield file.") ));
	return "";
}


// Custom Property - obtain an int
int Particle::getCustomProperty_Int(const std::string &_name) const {
	std::string data = getCustomProperty_String(_name);
	int value;
	if(str2int(data,value)!= 0 ){
		throw(ProcedureException( std::string("Atom ") + pdbname + 
			std::string(" ,property ") + 
			_name + 
			std::string(" is not an integer quantity.") ));
	}
	return value;
}

// Custom Property - obtain a double
double Particle::getCustomProperty_Double(const std::string &_name) const {
	std::string data = getCustomProperty_String(_name);
	double value;
	if(str2double(data,value)!= 0 ){
		throw(ProcedureException( std::string("Atom ") + pdbname + 
			std::string(" ,property ") + 
			_name + 
			std::string(" is not a floating point quantity.") ));
	}
	return value;
}





void Width4Printer( const std::string& _term )
{
	static StringBuilder sb;
	sb.setTo(_term);
	if(sb.size()<4)
	{
		sb.insert(0,' ');
	}
	sb.PadRight(4,' ');
	std::cout << sb;
}

void Particle::detail() const
{
	info( Verbosity::Loud, 4, true );
}

void Particle::info( Verbosity::Type verbosity, int padNumbers, bool endLine ) const
{
	//printf("ATOM   %4i  %-3s %3s   %3i    %8.3lf%8.3lf%8.3lf\n ",
	//	i, atom[i].pdbname.c_str(),
	//	atom[i].parentl3name.c_str(), atom[i].ir, atom[i].pos().x, atom[i].pos().y, atom[i].pos().z);

	if( verbosity > Verbosity::Silent )
	{
		Width4Printer( pdbname );
		std::cout << ' ' << std::setw(padNumbers) << i << ' ';
		Width4Printer( parentname );
		std::cout << ' ' << std::setw(padNumbers) << ir;
		if( endLine ) std::cout << std::endl;
		return; // just print enough info to identifty the particle
	}
	if( verbosity >= Verbosity::Normal )
	{
		Width4Printer( rawname );
		Width4Printer( pdbname );
		std::cout << ' ' << std::setw(padNumbers) << i << ' ';
		std::cout << ' ' << parentletter;
		Width4Printer( parentname );
		std::cout << ' ' << std::setw(padNumbers) << ir;
	}
	if( verbosity >= Verbosity::Loud )
	{
		std::cout  << ' ' << Z
			<< ' ' << mass
			<< ' ' << radius
			<< ' ' << epsilon
			<< ' ' << charge
			<< ' ' << epot;
	}
	//if( Verbosity::Type >= 3 )
	//{
	//	// Type/Status Info ??
	//}
	if( endLine ) std::cout << ' ' << std::endl;
}

size_t Particle::memuse(int level) const
{
	level--;
	return sizeof(Particle) +
		(cov12atom.size() +
		cov13atom.size() +
		cov14atom.size() )*sizeof(CovalentAtom);
}

void Particle::add12OnlyCovbond(int iatom2)
{
	// Prevent double addition by checking if this index is already present
	bool contains = false;
	for( size_t i = 0; i < cov12atom.size(); i++ )
	{
		if( cov12atom[i].i == iatom2 )
		{
			contains = true;
			break;
		}
	}
	if( !contains )
	{
		CovalentAtom covatom;
		covatom.i = iatom2;
		cov12atom.push_back(covatom);
	}
}

void Particle::offsetInternalIndices( int offset )
{
	for(unsigned i=0;i<cov12atom.size();i++)
	{
		cov12atom[i].i  += offset;
		cov12atom[i].i2 += offset;
		cov12atom[i].i3 += offset;
	}																			
	for(unsigned i=0;i<cov13atom.size();i++)
	{	
		cov13atom[i].i  += offset;
		cov13atom[i].i2 += offset;
		cov13atom[i].i3 += offset;
	}																			
	for(unsigned i=0;i<cov14atom.size();i++)
	{	
		cov14atom[i].i  += offset;
		cov14atom[i].i2 += offset;
		cov14atom[i].i3 += offset;
	}
}







// Residue member functions -------------------------------------------

size_t Residue::size() const
{
	ASSERT( ifirst >= 0 && ilast >= 0, UninitialisedException, "Residue object ifirst, ilast are not assigned");
	ASSERT( ilast >= ifirst, CodeException, "Residue object are miss-assigned: ilast < ifirst!");
	return ilast - ifirst + 1;
}

void Residue::addOffset(int offset)
{
	if( ifirst >= 0 ) ifirst += offset;
	if( ilast >= 0 ) ilast += offset;
	if( iCA >= 0 ) iCA += offset;
	if( iHA >= 0 ) iHA += offset;
	if( iN >= 0 ) iN += offset;
	if( iC >= 0 ) iC += offset;
	if( iO >= 0 ) iO += offset;
	if( iH >= 0 ) iH += offset;
	if( iCB >= 0 ) iCB += offset;
}

void Residue::info( Verbosity::Type verbosity ) const
{
	if( verbosity > Verbosity::Silent )
	{		
		if( param != NULL )
		{
			std::cout << param->letter;
			std::cout << ' ' << std::setw(4) << param->s_name();
		}	
		if(  verbosity >= Verbosity::Normal )
		{
			std::cout 
				<< "1st:" << std::setw(4) << ifirst << ' '
				<< "last:" << std::setw(4) << ilast << ' '
				<< "iCA:" << std::setw(4) << iCA << ' '
				<< "iHA:" << std::setw(4) << iHA << ' '
				<< "iN:" << std::setw(4) << iN << ' '
				<< "iC:" << std::setw(4) << iC << ' '
				<< "iO:" << std::setw(4) << iO << ' '
				<< "iH:" << std::setw(4) << iH << ' '
				<< "iCB:" << std::setw(4) << iCB;
		}
		std::cout << std::endl;
	}
}

