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

#ifndef __PROFILEBASE_H
#define __PROFILEBASE_H


namespace Protocol
{

	//-------------------------------------------------
	//
	/// \brief Base class for profiles for simple and more elaborate Temperature or Pressure
	/// Evolutions such as TempLinear, TempExp[onential], etc..
	///
	/// \details 
	///
	/// These classes provide a simple way to change the Temperature
	/// throughout a run of arbitrary nature. The Protocols can request
	/// the current Temperature by calling get(double progress ) where
	/// progress is a progress variable indicating how far through the run
	/// the current Temperature is.
	/// \author M.Tyka
	///
	///
	/// \author Mike Tyka 
	///
	///
	class PD_API ProfileBase{
	public:
		ProfileBase(){}
		virtual ~ProfileBase(){}
		virtual ProfileBase* clone() const = 0; 
		virtual double get(double progress) const = 0;
		virtual void info() const = 0;
	};






	//-------------------------------------------------
	//
	/// \brief Provides the most simple of Profiles, the constant profile. This could for example be a constant temperature or constant pressure 
	///        Not that this has nothing to do with thermo/pressure stating but instead simple provides a constant reference temperature or pressure. 
	/// \details 
	///    
	/// \author Mike Tyka 
	///
	class PD_API ConstantProfile: public ProfileBase {
	public:
		ConstantProfile(double _T0) { T0 = _T0; }
		virtual ~ConstantProfile(){}
		virtual ConstantProfile* clone() const 
		{ 
			return new ConstantProfile(*this); 
		}
		virtual double get(double progress) const {
			return T0;
		}
		virtual void info() const{
			printf("Constant Profile: %lf",T0);
		}
	protected:
		double T0;
	};






	//-------------------------------------------------
	//
	/// \brief Linear profile. I.e. a linear temperature or pressure ramp 
	/// \details 
	/// \author Mike Tyka 
	///
	class PD_API LinearProfile: public ConstantProfile {
	public:
		LinearProfile(double _T0, double _TN): ConstantProfile(_T0){
			TN = _TN;
		}
		virtual LinearProfile* clone() const 
		{ 
			return new LinearProfile(*this); 
		}
		virtual double get(double progress) const {
			return T0 + (progress * (TN - T0));
		}
		virtual void info() const{
			printf("Linear Profile: %lf --> %lf ",T0,TN);
		}
	protected:
		double TN;
	};







	//-------------------------------------------------
	//
	/// \brief Exponential profile. I.e. a exponential temperature or pressure ramp 
	/// \details 
	/// \author Mike Tyka 
	///
	class PD_API ExponentialProfile: public LinearProfile {
	public:
		ExponentialProfile(double _T0, double _TN): LinearProfile(_T0,_TN){
			alpha = -(log(T0 / TN));
		}
		virtual ExponentialProfile* clone() const 
		{ 
			return new ExponentialProfile(*this); 
		}
		virtual double get(double progress) const {
			return T0 * exp(alpha * progress);
		}
		virtual void info() const{
			printf("Exponential Profile: %lf --> %lf ",T0,TN);
		}
	protected:
		double alpha;
	};



}

#endif

