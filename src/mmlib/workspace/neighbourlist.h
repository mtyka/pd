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

#ifndef __NEIGHBOUR_LIST_H
#define __NEIGHBOUR_LIST_H

#include <valarray>
#include "componentbase.h"

#include "system/fundamentals.fwd.h"
#include "workspace/bondorder.fwd.h"
#include "workspace/workspace.fwd.h"

struct NeighbourData
{
	int n; // number of neighbors
	int *i; // index of neighbors
	int *Type; // Type
};







//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka  
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API NeighbourListBase: public WorkSpaceComponentBase
{
public:
	NeighbourListBase();
	virtual ~NeighbourListBase();

	virtual size_t  memuse(int level);

	// Neighbour Data Accessor
	inline const NeighbourData* getData() const { return fnbor; }

	// Configuration
	void     enable(){ Enabled = true; };
	void     disable(){ Enabled = false; };
	void     calcShadow( bool Enabled ) { CalcShadow = Enabled; }	
	void     requestCutoff( double newcutoff );
	void     setPadding( double newpadding );
	unsigned getFullUpdateCount() const { return m_FullUpdateCount; } 


	// OO Accessors
	virtual size_t   nNeighbours(size_t i) const  = 0;
	virtual int      getNeighbourIndex(size_t i, size_t nj) const  = 0;
	virtual int      getNeighbourBondOrder(size_t i, size_t nj) const  = 0;
	virtual int      getNeighbourImageNumber(size_t i, size_t nj) const  = 0;


	// Neighbor list calculations
	virtual void calcNewList() = 0; /// A full recalculation of the neighbour list

protected:
	// Internal memory management
	virtual void reinit( WorkSpace* _wspace );
	virtual void reassignList();
	virtual void reserveMemoryFor(int newmaxneighbors);
	void         clear();
	void         incFullUpdateCount() { m_FullUpdateCount++; }

	NeighbourListBase(const NeighbourListBase &copy);


protected: 	// protected member data ---------------------------------

	/// Should we perform the calculations at all?
	bool Enabled;  

	/// Cutoff to what get's included. 
	double Cutoff;  

	/// Padding to the Cutoff, i.e. atoms included in the neighbor list 
	double Padding;  

	/// included shadow neighbor in neighbor list
	bool CalcShadow;



	// Internal neighbourlist data allocation -----------------------------

	/// A large contiguos block of memory in which we store neighbour info.
	char *neighborlistspace;		

	/// A memory structure containing pointers to data within neighborlistspace
	NeighbourData *fnbor;				

	/// Indicates the 'currentMaxNeighbors' and hence the current amount of allocated memory to the pointers in fnbor
	int currentMaxNeighbors;

	/// This is a counter which gets updated on every full update of the neighbor 
	/// list. In this way, forcefields and other components can 
	/// tell if they must update themselves.
	unsigned m_FullUpdateCount;
};



// Accesory inline functions to decode the 32bit neighbor information:
inline unsigned int NList32Bit_Index(unsigned int rawdata){ return rawdata&0x00FFFFFF; };
inline unsigned int NList32Bit_BondOrder(unsigned int rawdata){ return (rawdata>>24)&7; };
inline unsigned int NList32Bit_ImageNumber(unsigned int rawdata){ return (rawdata>>27); };






//-------------------------------------------------
//
/// \brief  Base class for NeighborLists using 32Bits per neighbor 
///
///
/// \details The structure of the Neighbor int is:
/// // Bit structure:
/// 76543210 76543210 76543210 76543210
/// `---'`-' `-------< index >--------'
///   |    \------- bondorder            (0-8) 
///    \----------- image vector (index) (0-32) 
///
/// This type of NeighborList has a limit of 16.7 Million atoms
/// For larger systems the 64Bit Neighbour Lists must be used
///
/// \author Mike Tyka  
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API NeighbourList_32Bit_Base: public NeighbourListBase
{
public:
	NeighbourList_32Bit_Base():NeighbourListBase(){};
	virtual ~NeighbourList_32Bit_Base(){};
	virtual void reassignList();

	virtual size_t   nNeighbours(size_t i) const { return fnbor[i].n; }
	virtual int      getNeighbourIndex(size_t i, size_t nj) const {       return NList32Bit_Index(fnbor[i].i[nj]); }
	virtual int      getNeighbourBondOrder(size_t i, size_t nj) const {   return NList32Bit_BondOrder(fnbor[i].i[nj]); }
	virtual int      getNeighbourImageNumber(size_t i, size_t nj) const { return NList32Bit_ImageNumber(fnbor[i].i[nj]); }
	
protected:

	/// Reserves the memory to allow newmaxneighbors many 
	/// neighbors (in total, not per atom);
  virtual void reserveMemoryFor(int newmaxneighbors);

	/// Analyses the matrix of off-diagonal low bondorders and
  /// correct the bondorders of the respective neighbors. This is
  /// necessary for disulfide bonds for example
	virtual void correctOffDiagonalBondOrders();
};







//-------------------------------------------------
//
/// \brief  Basic 32Bit NeighborList for Infinite or PeriodicBox spaces
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka  
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API NeighbourList: public NeighbourList_32Bit_Base
{
public:
	NeighbourList():NeighbourList_32Bit_Base(){};
	virtual ~NeighbourList(){};
	// Neighbor list calculations
	virtual void calcNewList(); /// A full recalculation of the neighbour list

protected:
	virtual void calcNewList_InfiniteSpace();
	virtual void calcNewList_PeriodicBox();
};







//-------------------------------------------------
//
/// \brief  GroupBased 32Bit NeighborList for Infinite or PeriodicBox spaces
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka  
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API NeighbourList_GroupBased: public NeighbourList
{
public:
	NeighbourList_GroupBased():NeighbourList(){};
	virtual ~NeighbourList_GroupBased(){};
protected:
	virtual void calcNewList_PeriodicBox(); 
};




///// OLD DEPRECATED NEIGHBOR LISTS


/// NeighborList for InfiniteSpaces





//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka  
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API DeprecatedNeighbourList: public NeighbourListBase
{
public:
	DeprecatedNeighbourList():NeighbourListBase(){};
	virtual ~DeprecatedNeighbourList(){};


	virtual size_t   nNeighbours(size_t i) const { THROW(CodeException,"DEPRECATED CODE CALLED"); }
	virtual int      getNeighbourIndex(size_t i, size_t nj) const { THROW(CodeException,"DEPRECATED CODE CALLED"); }
	virtual int      getNeighbourBondOrder(size_t i, size_t nj) const { THROW(CodeException,"DEPRECATED CODE CALLED"); }
	virtual int      getNeighbourImageNumber(size_t i, size_t nj) const { THROW(CodeException,"DEPRECATED CODE CALLED"); }
	// Neighbor list calculations
	virtual void calcNewList(); /// A full recalculation of the neighbour list
};







//-------------------------------------------------
//
/// \brief  NeighborList for all types of spaces
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka  
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API NeighbourList_GeneralBoundary: public NeighbourListBase
{
public:
	NeighbourList_GeneralBoundary():NeighbourListBase(){};
	virtual ~NeighbourList_GeneralBoundary(){};


	virtual size_t   nNeighbours(size_t i) const { THROW(CodeException,"DEPRECATED CODE CALLED"); }
	virtual int      getNeighbourIndex(size_t i, size_t nj) const { THROW(CodeException,"DEPRECATED CODE CALLED"); }
	virtual int      getNeighbourBondOrder(size_t i, size_t nj) const { THROW(CodeException,"DEPRECATED CODE CALLED"); }
	virtual int      getNeighbourImageNumber(size_t i, size_t nj) const { THROW(CodeException,"DEPRECATED CODE CALLED"); }
	// Neighbor list calculations
	virtual void calcNewList(); /// A full recalculation of the neighbour list
};







//-------------------------------------------------
//
/// \brief  NeighborList for all types of spaces
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka  
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API NeighbourList_PeriodicBox: public NeighbourListBase
{
public:
	NeighbourList_PeriodicBox():NeighbourListBase(){};
	virtual ~NeighbourList_PeriodicBox(){};

	virtual size_t   nNeighbours(size_t i) const { THROW(CodeException,"DEPRECATED CODE CALLED"); }
	virtual int      getNeighbourIndex(size_t i, size_t nj) const { THROW(CodeException,"DEPRECATED CODE CALLED"); }
	virtual int      getNeighbourBondOrder(size_t i, size_t nj) const { THROW(CodeException,"DEPRECATED CODE CALLED"); }
	virtual int      getNeighbourImageNumber(size_t i, size_t nj) const { THROW(CodeException,"DEPRECATED CODE CALLED"); }


	// Neighbor list calculations
	virtual void calcNewList(); /// A full recalculation of the neighbour list
};





#endif

