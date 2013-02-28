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

#ifndef __BOUNDRY_H
#define __BOUNDRY_H

#include <valarray>

#include "system/system.h"





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
class PD_API Space 
{
public:
	Space(){};

	virtual void info() const = 0;

	virtual void getImage(
		Maths::dvector &imagedx,								///< image position
		int             ncell										///< cell number
	) const =0;
	
	virtual void getClosestImage(
		Maths::dvector &imagedx 								///< image position
	) const =0;

	virtual size_t ncells() const = 0;

	virtual void moveIntoBox(
		Maths::dvector &va
	) const = 0;

	virtual void moveIntoCell(
		Maths::dvector &va
	) const = 0;

	virtual void setupCellDimensions(double Cutoff) = 0;

	/// proportional scaling of box
	virtual void scale(double factor) = 0;

	virtual size_t         nBasisVectors() const { return 1; } 
	virtual Maths::dvector getBasisVector( size_t i) const { return Maths::dvector(0,0,0); } 
	virtual double         getSmallestChord() const { return -1; }

	double savedist;
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
class PD_API ClosedSpace: public Space 
{
public:
	ClosedSpace():Space(){};

	virtual double volume() const = 0;

	virtual void getEvenPointDistribution(std::vector<Maths::dvector> &pointarray, unsigned N) const = 0;

	// THIS SHOULD BE PURE VIRUTAL BUT FOR SOME REASON SWIG HIKCS UP. WILL SORT OUT LATER.

	virtual void getBoxVectors(Maths::dvector &A,Maths::dvector &B,Maths::dvector &C) const = 0; 
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
class PD_API InfiniteSpace: public Space 
{
public:
	InfiniteSpace(){};

	/// print information about the space
	virtual void info() const { printf("Space: InfiniteSpace\n"); }

	virtual void getImage(
		Maths::dvector &imagedx,								///< image position
		int             ncell										///< cell number
	) const {}
	
	virtual void getClosestImage(
		Maths::dvector &imagedx								///< image position
	) const {}

	virtual size_t ncells() const { return ncells_inline(); }

	inline size_t ncells_inline() const { return 1; }

	virtual void moveIntoBox(
		Maths::dvector &va
	) const  {}

	virtual void moveIntoCell(
		Maths::dvector &va
	) const  {}

	virtual void setupCellDimensions(double Cutoff){}

	virtual void scale(double factor){}

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
class PD_API PeriodicBox: public ClosedSpace 
{
public:
	/// creates a cubic PeriodicBox with specified dimension
	PeriodicBox(double allsize=1.0):  // space must have a default constructor
		boxSize(allsize,allsize,allsize)
	{
		init();
	};
	
	/// creates a PeriodicBox with specified dimensions
	PeriodicBox(const Maths::dvector &_boxSize): 
		boxSize(_boxSize)
		
	{	
		init();
	};

	/// creates a PeriodicBox with specified dimensions
	PeriodicBox(double x, double y,double z): 
		boxSize(Maths::dvector(x,y,z))
	{	
		init();
	};

	/// makes a box big enough to accomodate the system, adding on a margin
	PeriodicBox(System &sys, double margin, bool transformToOptimum=false);

	/// print information about the space
	virtual void info() const;

	// Generic virtual function - slow but general 
	virtual void getImage(
		Maths::dvector &imagedx,								///< image position
		int             ncell										///< cell number
	) const 
	{
		imagedx.add(celloffset[ncell]);
		while(imagedx.x >  halfCellSize.x) imagedx.x -= cellSize.x; 
		while(imagedx.x < -halfCellSize.x) imagedx.x += cellSize.x; 
		while(imagedx.y >  halfCellSize.y) imagedx.y -= cellSize.y; 
		while(imagedx.y < -halfCellSize.y) imagedx.y += cellSize.y; 
		while(imagedx.z >  halfCellSize.z) imagedx.z -= cellSize.z; 
		while(imagedx.z < -halfCellSize.z) imagedx.z += cellSize.z; 
	}

	// Generic virtual function - slow but general 
	virtual void getClosestImage(
		Maths::dvector &imagedx								///< image position
	) const 
	{
			while(imagedx.x >  halfBoxSize.x) imagedx.x -= boxSize.x; 
			while(imagedx.x < -halfBoxSize.x) imagedx.x += boxSize.x; 
			while(imagedx.y >  halfBoxSize.y) imagedx.y -= boxSize.y; 
			while(imagedx.y < -halfBoxSize.y) imagedx.y += boxSize.y; 
			while(imagedx.z >  halfBoxSize.z) imagedx.z -= boxSize.z; 
			while(imagedx.z < -halfBoxSize.z) imagedx.z += boxSize.z; 
	}


	virtual size_t ncells() const { return ncells_inline(); }
	inline size_t ncells_inline() const { return celloffset.size(); }

	virtual void moveIntoBox(
		Maths::dvector &va
	) const
	{
		getClosestImage(va);
	}

	virtual void moveIntoCell(
		Maths::dvector &va
	) const
	{
		while(va.x >  halfCellSize.x) va.x -= cellSize.x; 
		while(va.x < -halfCellSize.x) va.x += cellSize.x; 
		while(va.y >  halfCellSize.y) va.y -= cellSize.y; 
		while(va.y < -halfCellSize.y) va.y += cellSize.y; 
		while(va.z >  halfCellSize.z) va.z -= cellSize.z; 
		while(va.z < -halfCellSize.z) va.z += cellSize.z; 
	}

	virtual void setupCellDimensions(double Cutoff);


	virtual double volume() const 
	{
		return boxSize.x * boxSize.y * boxSize.z;
	}

	// provides a set of points even distributed through the space
	virtual void getEvenPointDistribution(std::vector<Maths::dvector> &pointarray, unsigned N) const;

	/// proportional scaling of box
	virtual void scale(double factor){
		boxSize.mul(factor);
		init();  // you must reinitialise !
	};

	virtual size_t nBasisVectors() const { return basisvector.size() ; }  
	virtual Maths::dvector getBasisVector(size_t i) const{ return basisvector[i];  } 
	virtual double getSmallestChord() const { return Maths::min( boxSize.x, boxSize.y, boxSize.z ); }

	virtual void getBoxVectors(Maths::dvector &A,Maths::dvector &B,Maths::dvector &C) const {
		A = Maths::dvector(boxSize.x,0,0);
		B = Maths::dvector(0,boxSize.y,0);
		C = Maths::dvector(0,0,boxSize.z);
	}

protected:
	void init(){
		halfBoxSize.setTo(boxSize);
		halfBoxSize.div(2.0);
		halfCellSize.setTo(boxSize);
		cellSize.setTo(boxSize);
		halfCellSize.div(2.0);

		vi.setTo(boxSize.x,0,0);
		vj.setTo(0,boxSize.y,0);
		vk.setTo(0,0,boxSize.z);

		savedist = Maths::min(boxSize.x/2,boxSize.y/2,boxSize.z/2);

		setupCellDimensions(1);
	}

public:
	/// vector from one corner to other corner through (0,0,0)
	Maths::dvector boxSize;
	Maths::dvector halfBoxSize;
	Maths::dvector cellSize;
	Maths::dvector halfCellSize;

	Maths::dvector vi,vj,vk;
	std::valarray<Maths::dvector> celloffset;
	std::valarray<Maths::dvector> basisvector;

	unsigned nx;
	unsigned ny;
	unsigned nz;
	
	double m_cutoff;
	double m_sqrcutoff;
	bool   need_replicas;
};

// converters


// just does a dynamic_cast and retruns the result. I.e. if the Space passed to 
// this function does not derive from ClosedSpace the result will be NULL!
ClosedSpace *castTo_ClosedSpace( Space *sp );


#endif


