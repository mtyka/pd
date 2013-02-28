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

#ifndef __MONEXP_H
#define __MONEXP_H

#include "monitors/monitorbase.h" // Provides a base class
#include "forcefields/forcefield.h"
#include "workspace/workspace.fwd.h"

namespace Monitors{


	//-------------------------------------------------
	//
	/// \brief Monitors Energy Differences calculated elsewhere and
	/// takes the FEP (Free Energy MoveBase) exponential as it's data
	///
	/// \author Mike Tyka 
	///

	class PD_API FEP_Monitor: public Monitors::MonitorBase{
	public:
		FEP_Monitor(): MonitorBase() {
			name = "FEP";
			Temperature = 300;
		}

		FEP_Monitor( double * new_epot_diff): MonitorBase() {
			name = "FEP";
			Temperature = 300;
			add_epot_diff(new_epot_diff);
		}

		~FEP_Monitor(){

		};

		virtual void setcurdata(){
			double total_energy_difference = 0.0;

			for(size_t i=0;i < energy_difference.size(); i++ ){
				total_energy_difference +=
					*(energy_difference[i]);
			}

			addData(  exp(-total_energy_difference/(Physics::PhysicsConst::kB*Temperature)) );
		}

		void add_epot_diff(double * new_epot_diff){
			energy_difference.push_back( new_epot_diff );
		}


		// Modified Averaging function that return -RT ln( AverageData ) instead
		// of just data

		double fe_average()const {
			return fe_average(0,nData());
		}
		double fe_average(size_t start)const {
			return fe_average(start,nData());
		}

		double fe_average(size_t start, size_t end)const {
			return -Temperature*Physics::PhysicsConst::R_kcal * log( getAverageSubdata(start,end) );
		}


		double fe_average_prop(double start, double end)const {
			if((start<0.0)||(start>1.0)||
				(end<0.0)||(end>1.0)){
					THROW(ArgumentException,"ERROR: start and end must be between 0.0 and 1.0");
			}

			return fe_average(unsigned(start*double(nData())),
				unsigned(end*double(nData())));
		}

		/// Performs a bootstrap block analysis of the free energy averages etc..
		void fe_bootstrap(double ignore, unsigned blocks, int blocks_end=-1);

		/// Performs a bootstrap block analysis of the free energy averages and returns standard deviation
		double fe_bootstrap_sd(double ignore, unsigned nblocks);

		void fe_print_running_av(unsigned start, unsigned end, unsigned Step=1)const {

			for(unsigned finish = (start+Step); finish < end; finish += Step ){
				printf("%s %5d %5d %.7e\n",name.c_str(), start, finish, fe_average(start,finish) );
			}
		}


		void fe_print_running_av_prop(double start, double end, unsigned Step=1)const {
			if((start<0.0)||(start>1.0)||
				(end<0.0)||(end>1.0)){
					THROW(ArgumentException,"ERROR: start and end must be between 0.0 and 1.0");
			}

			return fe_print_running_av(unsigned(start*double(nData())),
				unsigned(end*double(nData())),
				Step);
		}

		double Temperature;
	protected:
		std::vector< double* > energy_difference;

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

	class PD_API FEP_ResConf_Cartesian:public Monitors::MonitorBase{
	public:

		FEP_ResConf_Cartesian(WorkSpace *_wspace): MonitorBase() {
			settodefault();
			wspace = _wspace;
			eqpos = NULL;
		}
		~FEP_ResConf_Cartesian(){
			delete[]eqpos;
		};

		double k;
		double Power;
		int Type;
		bool addtototalene;
		int FirstRes;
		int LastRes;

		virtual void settodefault(){
			k = 0;
			Power = 2;
			Type = 0;
			addtototalene = true;
			FirstRes =-1;
			LastRes = -1;
			ncalls = 0;
			nskip=0;
		};

		virtual void setcurdata();

		virtual double calcEnergyAtQ(double newQ);

		virtual void info() const; // prints a little block of parameter information

		void saveCurrentAtomPositions();

		virtual void printAllData() const;
		virtual int setup();

		int nskip;
	protected:
		virtual void infoLine() const; // prints a line of current energies
		virtual void infoLineHeader() const; // prints the headers for the above function
		Maths::dvector *eqpos; // saves the equilibrium atom positions;

		std::vector <double> res_a_fep;
		std::vector < std::vector<double> > res_ab_fep;
		std::vector < std::vector< std::vector<double> > > res_abc_fep;

		int ncalls;
		int ndatapoints;
	private:
		WorkSpace *wspace;
	};


}

#endif
