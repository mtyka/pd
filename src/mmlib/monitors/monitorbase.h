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

#ifndef __MONITORBASE_H
#define __MONITORBASE_H

// Essential Headers
#include "object.h" // Provides a base class
#include "maths/histogram.h" // Provides a class member
#include "workspace/workspace.fwd.h"

// Forward declarations
namespace Physics
{
	class PD_API ForcefieldBase;
	class PD_API RestraintForcefieldBase;
	class PD_API FF_Extension_TI;
}

class PD_API MoleculeBase;

namespace Monitors
{
	/// \brief Abstract base class for all monitors. Implements a number of generic
	/// \brief Functions that allow dealing with data output & averaging
	/// \author M.Tyka
	///
	/// A Monitor is a class that performs some sort of measurement on a WorkSpace,
	/// Forcefield or any other property of a simulation system. Usually measurements
	/// need to taken at regular intervals through a simulation in order to obtain
	/// ensemble averages although Monitors can also be used for single measurements.
	/// The Interface is layed out in the abstract class MonitorBase and defines an
	/// abstract function setcurdata() which derived classes must overload, implementing
	/// their particular measurements. Monitors can be added to a simulation protocol
	/// by calling Protocol::addMonitor() which will add a pointer to the respective
	/// Monitor into an internal MonitorContainer class, ”registering” it for execution
	/// during the protocol. When run, the protocol will call the measure() function of
	/// every stored monitor measure() takes the data calculated by the overloaded
	/// function setcurdata() and stores it away in a data array. Following the
	/// simulation, MonitorBase provides a large variety of function that allow
	/// startistical analysis of the stored data such as calculation of averages,
	/// generation of histograms and bootstrap block average analyses. This object
	/// oriented design makes implementation of new monitors trivial, since only a
	/// single function has to be overloaded while the full breadth of functionality
	/// is inherited from MonitorBase.
	///
	/// Monitors implement all kinds of Passive measurements done during a or multiple
	/// simulations. The abstract base class MonitorBase implements a large number of useful functions
	/// for collection, storing and analysing the measured data.
	/// Implemented Classes derived from MonitorBase must implement a function named
	/// setcurdata() - this function is called automatically by the base class function
	/// measure() every time a measurement is to occur. measure() is usually called by a Class
	/// derived from ProtocolBase. For example a class PD_API implementing a Molecular Dynamics
	/// Protocol could call measure() after every time Step.
	/// measure() will then call setcurdata() which must perform the measurement and
	/// call addData(double data) to add it to the stack or store of data.
	/// into a std::vector array for storage.
	///
	/// At the end of the run MonitorBase provides a number of helper functions to
	/// make statistical analyses of the data, create histograms, print & write the data.
	/// External functions provide ways to combine data from multiple monitors to create 2D
	/// histograms or WHAM analyses.

	class PD_API MonitorBase : public Object
	{
	public:
		MonitorBase();
		virtual MonitorBase* clone() const = 0;

		virtual int measure();
		virtual int repeat();
		virtual void reset();
		virtual void printHeader() const;
		virtual void printCurData() const;
		virtual void printData(size_t index) const;
		virtual void printAllData() const;
		virtual void printHistogram(double binwidth);
		virtual void printHistogram(double binwidth, unsigned start, unsigned end) const;

		/// returns the number of data elements currently stored up
		size_t nData() const{ return dataarray.size(); }

		/// returns last piece of data
		double getCurData() const { return curdata; };

		/// adds a piece of data to the data stack
		virtual int  addData(double data);

		/// returns a particular piece of data with index 'index'
		double getData(size_t index) const;

		/// returns arithmetic average of data between indices start and (including) end
		virtual double getAverageSubdata(size_t start, size_t end) const;

		/// returns arithmetic average of all data
		virtual double getAverageSubdata() const;

		/// returns arithmetic average from 'start' to end of data
		virtual double getAverageSubdata(size_t start) const;

		/// returns variance of data in range
		virtual double getVarianceSubdata(size_t start, size_t end) const;

		///\brief returns arithmetic average of proportion of data bounded by start and double
		///
		/// Both parameters must be between 0.0 and 1.0
		virtual double getVariance(double start=0.0, double end=1.0) const;

		///\brief returns arithmetic average of proportion of data bounded by start and double
		///
		/// Both parameters must be between 0.0 and 1.0
		virtual double getAverage(double start=0.0, double end=1.0) const;

		/// Performs a Bootstrap block average analysis
		virtual void printBlockAverage(double ignore, unsigned blocks, int blocks_end=-1);

		/// Performs a Bootstrap block average analysis, returns standard deviation
		virtual double getBlockAvStdDev(double ignore, unsigned blocks);

		/// prints a running average from start to end
		virtual void printRunningAverageSubdata(unsigned start, unsigned end, unsigned Step=1) const;

		/// prints a running average from proportion start to end
		virtual void printRunningAverage(double start=0.0, double end=1.0, unsigned Step=1) const;

		/// Returns a 1D histogram of the data with a given binwidth , ignoring data before index 'start'
		Maths::Histogram1D getHistogram(double binwidth, unsigned start=0) const;

		/// Returns a 1D histogram of the data within a given range
		Maths::Histogram1D getHistogram(double binwidth,
			unsigned start,
			unsigned end) const;

		// Reweighted functions

		/// returns arithmetic average of data reweighted by weights from 'weight'
		virtual double getAverageSubdata(const MonitorBase &weight, unsigned start, unsigned end) const;

		/// returns reweighted arithmetic average of all data
		virtual double getAverageSubdata(const MonitorBase &weight) const;

		/// returns reweighted arithmetic average from 'start' to end of data
		virtual double getAverageSubdata(const MonitorBase &weight,unsigned start) const;

		///\brief returns reweighted arithmetic average of proportion of data bounded by start and double
		///
		/// start&end parameters must be between 0.0 and 1.0
		virtual double getAverage(const MonitorBase &weight, double start=0.0, double end=1.0) const ;

		/// Performs a bootstrap analysis on the reweighted data
		virtual void printBlockAverage(const MonitorBase &weight,double ignore, unsigned blocks, int blocks_end=-1);

		/// Performs a bootstrap analysis on the reweighted data & returns standard deviation
		virtual double getBlockAvStdDev(const MonitorBase &weight,double ignore, unsigned nblocks);

		/// Returns a histogram of the reweighted data.
		Maths::Histogram1D_double getHistogram(const MonitorBase &weight,double binwidth, unsigned start=0) const;
		Maths::Histogram1D_double getHistogram(const MonitorBase &weight,
			double binwidth,
			unsigned start,
			unsigned end) const;

		// public data (parameters)
		Verbosity::Type OutputLevel;									///< completely !OutputLevel mode (no screen output at all)
		bool PrintAverage;						///< also print running average

		unsigned RawIgnoreFirstN;
		unsigned RawSkip;

		void setGranularity(unsigned newgran);

	protected:
		virtual void setcurdata()= 0;		///< derived classes must override this with their particular calculation

		double curdata;           ///< The most recent piece of data
		std::string formatstring;				///< This can be overriden by derived classes if they desire a different output format
		std::string headerformatstring; ///< This can be overriden by derived classes if they desire a different output format

	private:

		std::vector< double > dataarray; ///< Contains the entire data of the run
		unsigned granularity;
		double gran_buffer;
		double gran_count;
		double lastdata;
		bool   lastdata_valid;
	protected:
		unsigned raw_count;
	};
}

#ifdef SWIG
%template(ObjectContainer_Monitor) ObjectContainer<Monitors::MonitorBase>;
#endif

namespace Monitors
{

	class PD_API MonitorContainer: public MonitorBase, public ObjectContainer<MonitorBase>
	{
	public:
		MonitorContainer();
		virtual MonitorContainer* clone() const { return new MonitorContainer(*this); }

		virtual int measure();
		virtual int repeat();
		virtual void reset();
		virtual void printHeader() const;
		virtual void printCurData() const;
		virtual void printAllData() const;
		virtual void printAllDataParallel() const;

	private:
		/// A container doesnt actually make any measurements of it's own so
		/// this function is technically irelevant in this class. Thus measure()
		/// is called which then calls measure() in all the child monitors.
		/// However MonitorContainer::measure() should really be called directly
		/// so this function here is private so that users dont call it directly.
		/// (it is nessessary tho cos MonitorContainers can contain Containers themselves)
		virtual void setcurdata()
		{
			measure();
		}
	};
}

#endif
