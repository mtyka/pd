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

#ifndef __OBJECT_H
#define __OBJECT_H

#include <string>
#include <vector>
#include "interface.h"

/// Object is a low level abstract base class other things can derive from
/// to inherit properties Object provides: The ability to compare the
/// identities of two objects in a better way than comparing their pointer addresses;

class PD_API Object : public ICloneable
{
public:	

	/// Default Constructor	
	Object(); 

	/// Copy constructor
	Object( const Object & _clone ); 
	virtual ~Object(){}
	
#ifndef SWIG
	/// Assignement operator. Instead of making an exact clone, we assign a unique
	/// object_id as if making anew object.
	Object &operator= ( const Object & _clone );
	friend bool operator== (const Object &a, const Object &b);
#endif


	/// generic identifier to be freely used by user;
	std::string name; 

	long getid() const { return object_id; }

private:

	/// Derived classes should override this
	virtual void setInternalName(); 

	/// static variable counting up unique object identifiers
	static long object_nextid;  

	/// unique id of this object
	long object_id;             
};

#ifndef SWIG
bool operator== (const Object &a, const Object &b);
#endif

/// A special Type of container:
///
/// Unlike the STL containers this contaner can store
/// Polymorphic Objects derived from a common baseclass PD_API T in
/// two distinct ways:
/// a) as copied, owned object in the same way std::vector does -
/// the functions addWithOwnership() and push_back() are used to add (or copy) an object
/// to the end of the list - a new instance is made on addition whose
/// ownership remains with this class. The instance is deleted when the
/// container dies. For this Type of Insertion the object class PD_API must provide
/// a function clone() which must return a pointer to a copied instance of itself.
///
/// b) as lent out pointers, owned by the caller. These can be inserted
/// using add() passing a pointer to a user-owned object. ObjectContainer
/// will not delete the object the pointer points to - that is the responsibility
/// of the caller. It is also the responsability of the user to remove the object
/// from the list before it is deleted! This can be done using erase()
///
/// NOTE: All stored classes must support ICloneable

template <class T>
class ObjectContainer
{
public:
	ObjectContainer()
	{
		// these flags can be used by derived classes to monitor if the contents
		// of their elements have changed and react accordingly
		added_element = false;
		removed_element = false;
		changed_element = false;
	}

	virtual ObjectContainer* clone() const
	{
		return new ObjectContainer(*this);
	}

	virtual ~ObjectContainer()
	{
		clear();
	}

	void clear()
	{
		removed_element = true; // flag that the elements have changed
		if( data.size() == 0 ) return;
		size_t datasize = data.size();
		for(size_t i=0;i<datasize;i++) pop_back();
		data.clear();
		data_owned.clear();		
	}

	ObjectContainer( const ObjectContainer& clone )
	{
		*this = clone;
	}

#ifndef SWIG
	ObjectContainer& operator=( const ObjectContainer& clone )
	{
		ASSERT( clone.data.size() == clone.data_owned.size(), CodeException, "Internal copy error");

		data.clear();
		data_owned.clear();

		for( size_t i = 0; i < clone.data.size(); i++ )
		{
			if( clone.data_owned[i] )
			{
				// Copy any data that the clone owns
				data.push_back( clone.data[i]->clone() );
				data_owned.push_back( true );
			}
			else
			{
				// We assume that the object is a 'create once' and that it is ok to pass round pointers to it.
				// These are usually python objects or objects created at the begining of program execution.
				data.push_back( clone.data[i] );
				data_owned.push_back( false );
			}
		}

		return *this;
	}
#endif 

	size_t size() const { return data.size(); }

	/// \brief Add an element to the container. 
	/// Note that ownership of the object
	/// remains with the caller who is responsible for 
	/// a) mainainting the object alive and in scope for the lifetime of the 
	///    container
	/// b) deleting the object *after* the container goes out of scope.

	void add( T & newelement )
	{
		preAdditionHook(newelement);
		data.push_back( &newelement );
		data_owned.push_back(false);
		added_element = true; // flag that the elements have changed
		postAdditionHook(newelement);
	}

protected:
	/// This function is empty, but can be overloaded by deriving classes
	/// to implement some kind of argument checking. In case of an error,
	/// an exception can be thrown
	virtual void preAdditionHook( T & newelement ){};

	/// This function is empty, but can be overloaded by deriving classes
	/// to implement some kind of post-argument-addition code. 
	virtual void postAdditionHook( T & newelement ){};

public:

	// Do not export this to python, as it's very incompatible with swig's
	// internal garbage collection. You should always use add() from python
	// in which case python retains the owner ship. 
#ifndef SWIG
	/// Add an element but transfer the owner ship to the container, which 
	/// will delete the object at the end of it's life time.
	void addWithOwnership( T* newelement )
	{
		ASSERT( newelement != NULL, NullArgumentException, "addWithOwnership() must not be given NULL elements");
		preAdditionHook(*newelement);
		data.push_back( newelement );
		data_owned.push_back(true);
		added_element = true; // flag that the elements have changed
		postAdditionHook(*newelement);
	}
#endif

	void erase(size_t index)
	{
		if((index>=data.size())) THROW (OutOfRangeException, "ObjectCOntainer: Index is out of bounds on erase() call");
		if( data_owned[index] )
		{
			delete data[index];
		}
		{
			typename std::vector< T* >::iterator iter = data.begin();
			iter += index;
			data.erase(iter);
		}
		{
			std::vector< bool >::iterator iter = data_owned.begin();
			iter += index;
			data_owned.erase(iter);
		}
		removed_element = true; // flag that the elements have changed
	}

	void pop_back()
	{
		ASSERT( data.size() > 0, OutOfRangeException, "ObjectCollection: Array must contain data for pop_back to be called.");
		size_t indexer = data.size()-1;
		if( data_owned[indexer] )
		{
			T* delData = data[indexer];
			delete delData;
		}
		data.pop_back();
		data_owned.pop_back();
		removed_element = true; // flag that the elements have changed
	}

#ifndef SWIG
	T &operator[] (size_t _element)
	{
		D_ASSERT( _element < data.size(), OutOfRangeException, "ObjectCOntainer: Index is out of bounds on operator[] call");
		return *(data[_element]);
	}
#endif

	T &element(size_t _element)
	{
		D_ASSERT( _element < data.size(), OutOfRangeException, "ObjectContainer: Index is out of bounds on element() call");
		return *(data[_element]);
	}

	const T &element(size_t _element) const
	{
		D_ASSERT( _element < data.size(), OutOfRangeException, "ObjectContainer: Index is out of bounds on const element() call");
		return *(data[_element]);
	}

protected:

	/// set by this base class
	bool added_element; 

	/// set by this base class
	bool removed_element; 

	/// not set by base class, only provided
	bool changed_element; 

private:

	/// pointers to the data
	std::vector< T* > data; 

	/// do we own the data ?
	std::vector< bool > data_owned; 
};

#endif




