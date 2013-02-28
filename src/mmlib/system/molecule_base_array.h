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


#ifndef __MOLECULE_BASE_ARRAY_H
#define __MOLECULE_BASE_ARRAY_H

// System Includes
#include <stdlib.h>
#include <string>
#include <vector>

//-------------------------------------------------
//
/// \brief  This class encapsulates the uber-secret-mega-private data of the  
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

template <class T>
class MoleculeBase_Array_Data
{
protected:
	MoleculeBase_Array_Data(){} // protected constructur prevents instantiation

	// ---[ super private data ]------------------------------------
	T * m_Data;			// stores the data
	size_t       m_N_Data;		// how many valid m_Datas
	size_t       m_N_Allocated;	// how many m_Datas is there space for
};


//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details 
/// This is a specialized private class for usage only by Molecule,
/// and Workspace. It's interface is based on that of std::vector but the 
/// implementation is much simpler and doesnt support special features such as
/// allocators etc.. 
///
/// \author Mike Tyka
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
template <class T>
class MoleculeBase_Array: public MoleculeBase_Array_Data< T > 
{
public:
	friend class PD_API MoleculeBase;
	friend class PD_API Molecule;
	friend class PD_API WorkSpace;

private:
	// constructors - means only friend can actually have one of these special
	// arrays, as intended by design (it would be useless to anyone else anyway)

	// Set up empty
	MoleculeBase_Array();

	// Destruct & clean memory
	~MoleculeBase_Array();

	/// 
	MoleculeBase_Array(const MoleculeBase_Array & da);

	/// Copy constructor 
	const MoleculeBase_Array<T> &operator=(const MoleculeBase_Array<T> & da);

public:
	// These functinos are publically available to allow anyone to inspect
	// the contents of this container
#ifndef SWIG
	/// Returns constant reference to an m_Data 
	const T & operator[](size_t index) const
	{
		D_ASSERT( index < MoleculeBase_Array_Data<T>::m_N_Data, OutOfRangeException, "" );
		return MoleculeBase_Array_Data<T>::m_Data[index];
	}
#endif
	/// Returns constant reference to an m_Data 
	const T & at(size_t index) const
	{
		D_ASSERT( index < MoleculeBase_Array_Data<T>::m_N_Data, OutOfRangeException, "" );
		return MoleculeBase_Array_Data<T>::m_Data[index];
	}

	size_t size() const 
	{
		return MoleculeBase_Array_Data<T>::m_N_Data;
	}

	//private:
#ifndef SWIG
	T & operator[](size_t index) 
	{
		D_ASSERT( index < MoleculeBase_Array_Data<T>::m_N_Data, OutOfRangeException, "" );
		return MoleculeBase_Array_Data<T>::m_Data[index];
	}
#endif

	T & at(size_t index) 
	{
		D_ASSERT( index < MoleculeBase_Array_Data<T>::m_N_Data, OutOfRangeException, "" );
		return MoleculeBase_Array_Data<T>::m_Data[index];
	}

	// ---[ Functions that change the size of the array ]-------
	// The following are functions used by the friends of this container
	// to modify/edit etc contents. These friends may or may not supply
	// wrappers to allow their users to do the editing but that's up to 
	// them 

private:
	/// Clears all m_Datas in this array
	void clear() {MoleculeBase_Array_Data<T>::m_N_Data = 0;};

	/// trims the memory allocation to the minimum 
	void trim();	

	/// adds an m_Data at the end (using references)
	int  push_back(const T & _newData);	

	/// removes the last m_Data
	void pop_back();		                  

	/// removes an arbitary m_Data
	void erase(size_t index);
};

template < class T >
MoleculeBase_Array < T >::MoleculeBase_Array()
{
	MoleculeBase_Array_Data<T>::m_N_Allocated = 1;
	MoleculeBase_Array_Data<T>::m_Data = new T[MoleculeBase_Array_Data<T>::m_N_Allocated];
	MoleculeBase_Array_Data<T>::m_N_Data = 0;
}

template < class T >
MoleculeBase_Array < T >::~MoleculeBase_Array()
{
	delete[] MoleculeBase_Array_Data<T>::m_Data;
}

template < class T >
MoleculeBase_Array < T >::MoleculeBase_Array (const MoleculeBase_Array < T > & da)
{
	MoleculeBase_Array_Data<T>::m_Data = NULL;
	(*this) = da;
}

template < class T >
const MoleculeBase_Array<T>  &MoleculeBase_Array < T >::operator=(const MoleculeBase_Array < T > & da)
{
	MoleculeBase_Array_Data<T>::m_N_Allocated = da.MoleculeBase_Array_Data<T>::m_N_Allocated;
	delete [] MoleculeBase_Array_Data<T>::m_Data;
	MoleculeBase_Array_Data<T>::m_Data = new T[MoleculeBase_Array_Data<T>::m_N_Allocated];
	MoleculeBase_Array_Data<T>::m_N_Data = da.MoleculeBase_Array_Data<T>::m_N_Data;
	for(size_t i = 0; i < MoleculeBase_Array_Data<T>::m_N_Data; i++) 
	{	
		// copy data using copy constructor
		MoleculeBase_Array_Data<T>::m_Data[i] = da.MoleculeBase_Array_Data<T>::m_Data[i];
	}
	return (*this);
}

// cuts the excess memory usage
template < class T >
void MoleculeBase_Array < T >::trim()
{
	// make new memory
	T *tempMemory = NULL;
	// minimise allocated memory 
	MoleculeBase_Array_Data<T>::m_N_Allocated = MoleculeBase_Array_Data<T>::m_N_Data;
	tempMemory = new T[MoleculeBase_Array_Data<T>::m_N_Allocated];
	if(tempMemory == NULL) THROW( OutOfMemoryException, "Out of memory while allocating in MoleculeBase_Array" );
	// transfer data
	for(size_t i = 0; i < MoleculeBase_Array_Data<T>::m_N_Data; i++) 
	{	
		// copy data using copy constructor
		tempMemory[i] = MoleculeBase_Array_Data<T>::m_Data[i];
	}
	delete[]MoleculeBase_Array_Data<T>::m_Data;
	MoleculeBase_Array_Data<T>::m_Data = tempMemory;
}

template < class T >
int MoleculeBase_Array < T >::push_back(const T & _newData)
{
	if(MoleculeBase_Array_Data<T>::m_N_Data >= 
		MoleculeBase_Array_Data<T>::m_N_Allocated) 
	{	
		// run out of space, change
		// make new memory
		T *tempMemory = NULL;
		if( MoleculeBase_Array_Data<T>::m_N_Allocated < 6){
			MoleculeBase_Array_Data<T>::m_N_Allocated += 1;
		}else{
			MoleculeBase_Array_Data<T>::m_N_Allocated *= 2;
		}
		tempMemory = new T[MoleculeBase_Array_Data<T>::m_N_Allocated];
		if(tempMemory == NULL) THROW( OutOfMemoryException, "Out of memory while allocating in MoleculeBase_Array" );
		// transfer data
		for(size_t i = 0; i < MoleculeBase_Array_Data<T>::m_N_Data; i++) 
		{	
			// copy data using copy constructor
			tempMemory[i] = MoleculeBase_Array_Data<T>::m_Data[i];
		}
		delete[]MoleculeBase_Array_Data<T>::m_Data;
		MoleculeBase_Array_Data<T>::m_Data = tempMemory;
	}

	MoleculeBase_Array_Data<T>::m_Data[MoleculeBase_Array_Data<T>::m_N_Data] = _newData;
	MoleculeBase_Array_Data<T>::m_N_Data++;

	return 0;
}

template < class T >
void MoleculeBase_Array < T >::pop_back()
{
	MoleculeBase_Array_Data<T>::m_N_Data--;
}

template < class T >
void MoleculeBase_Array < T >::erase(size_t index)
{
	for(int i = index; i < (MoleculeBase_Array_Data<T>::m_N_Data - 1); i++) 
	{	
		// move data using copy constructor
		MoleculeBase_Array_Data<T>::m_Data[i] = MoleculeBase_Array_Data<T>::m_Data[i + 1];
	}
	MoleculeBase_Array_Data<T>::m_N_Data--;
}

typedef MoleculeBase_Array<Particle> ParticleStore;
typedef MoleculeBase_Array<Residue> ResidueStore;
typedef MoleculeBase_Array<MoleculeRange> MoleculeRangeStore;


#endif


