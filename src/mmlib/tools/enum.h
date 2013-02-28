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

#ifndef __ENUM_TOOLS
#define __ENUM_TOOLS







//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details 
/// Code taken from:
/// http://www.codeproject.com/cpp/InheritEnum.asp
/// By Jon Rea of use in the MMLib project.
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
template <typename EnumT, typename BaseEnumT>
class InheritEnum
{
public:
  InheritEnum() {}
  InheritEnum(EnumT e)
	: enum_(e)
  {}

  InheritEnum(BaseEnumT e)
	: baseEnum_(e)
  {}

  explicit InheritEnum( int val )
	: enum_(static_cast<EnumT>(val))
  {}

  operator EnumT() const { return enum_; }
private:
  // Note - the value is declared as a union mainly for as a debugging aid. If
  // the union is undesired and you have other methods of debugging, change it
  // to either of EnumT and do a cast for the constructor that accepts BaseEnumT.
  union
  {
	EnumT enum_;
	BaseEnumT baseEnum_;
  };
};

#endif

