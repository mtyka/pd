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

#include "monexp.h"

#include "workspace/workspace.h"

using namespace Maths;
using namespace Physics;

namespace Monitors{

	/// Performs a bootstrap block analysis of the averages etc..
	void FEP_Monitor::fe_bootstrap(double ignore,
		unsigned blocks,
		int blocks_end){
			unsigned nblocks;

			argcheck_range("ignore",ignore,0.0,1.0);
			argcheck_gt("blocks",blocks,(unsigned)1);

			if(((int)blocks_end) < ((int)blocks)) blocks_end = blocks;
			if(blocks_end == blocks){
				printf("bs %s av % 8.4lf | ",name.c_str(), fe_average_prop(ignore,1.0) );
			}else{
				printf("bs %s av % 8.4lf\n",name.c_str(), fe_average_prop(ignore,1.0) );
			}

			for(nblocks = blocks; nblocks <= blocks_end; nblocks++){ // for various number of blocks
				if(blocks_end != blocks){
					printf("bs %s %2d ",name.c_str(), nblocks);
				}

				double sumx=0.0;
				double sumx2=0.0;
				unsigned b;
				for(b=0;b<nblocks; b++){ // calculate average over each block
					double bav = fe_average_prop(
						Maths::min(1.0,ignore + double(b)*(1.0 - ignore)/double(nblocks) ) ,
						Maths::min(1.0,ignore + double(b+1)*(1.0 - ignore)/double(nblocks) ) );
					printf("% 8.4lf ", bav );
					sumx += bav;
					sumx2 += Maths::sqr(bav);
				}
				for(;b<blocks_end;b++){
					printf(" ");
				}
				printf("stddev: % 8.4lf", sqrt( (sumx2 - Maths::sqr(sumx)/double(nblocks))/(nblocks-1) ) );
				printf("\n");
			}

	}

	/// Performs a bootstrap block analysis of the averages etc..
	double FEP_Monitor::fe_bootstrap_sd(double ignore, unsigned nblocks){
		argcheck_range("ignore",ignore,0.0,1.0);
		argcheck_gt("blocks",nblocks,(unsigned)1);

		double sumx=0.0;
		double sumx2=0.0;
		unsigned b;
		for(b=0;b<nblocks; b++){ // calculate average over each block
			double bav = fe_average_prop(
				Maths::min(1.0,ignore + double(b)*(1.0 - ignore)/double(nblocks) ) ,
				Maths::min(1.0,ignore + double(b+1)*(1.0 - ignore)/double(nblocks) ) );
			sumx += bav;
			sumx2 += Maths::sqr(bav);
		}
		return sqrt( (sumx2 - Maths::sqr(sumx)/double(nblocks))/(nblocks-1) );
	}

	void FEP_ResConf_Cartesian::info() const{ // prints a little block of parameter information
		printf(" --FEP_ResConf_Cartesian------------ \n");
		printf(" k: %8.4lf kcal/mol/A^2\n", k);
		printf(" Type: %4d \n", Type);

	}

	void FEP_ResConf_Cartesian::infoLine() const {
		printf("\t-----");
	}

	void FEP_ResConf_Cartesian::infoLineHeader() const {
		printf("\t-----");
	}

	int FEP_ResConf_Cartesian::setup(){

		delete [] eqpos;
		eqpos = new Maths::dvector[wspace->atom.size()];

		// This is already done by saveCurrentAtomPositions();
		//for(int i = 0; i < wspace->atom.size(); i++)
		//	eqpos[i].setTo(0, 0, 0, -1);

		saveCurrentAtomPositions();

		res_a_fep.clear();
		for(size_t ir=0;ir< wspace->res.size();ir++){
			res_a_fep.push_back(0.0);
		}
		res_ab_fep.clear();
		for(size_t jr=0;jr< wspace->res.size();jr++){
			res_ab_fep.push_back(res_a_fep);
		}
		res_abc_fep.clear();
		for(size_t jr=0;jr< wspace->res.size();jr++){
			res_abc_fep.push_back(res_ab_fep);
		}
		ndatapoints = 0;
		return 0;
	}

	void FEP_ResConf_Cartesian::saveCurrentAtomPositions(){
		printf("Harmonic Contraint: Reloading atom positions \n");
		for(int i = 0; i < wspace->atom.size(); i++) {
			eqpos[i].setTo(wspace->cur.atom[i].p);
		}
	}

	void FEP_ResConf_Cartesian::setcurdata(){
		int i;
		double ene;
		dvector force;

		double k_SI = k / (PhysicsConst::J2kcal * PhysicsConst::Na); // convert to SI units for calculation
		double forcemag;
		double dist;

		if(eqpos == NULL) {
			printf("FORCEFIELD ERROR: Restraint forcefield: Equilibrium positions have not been assigned \n");
			return;
		}

		std::vector<double> resene;
		ncalls++;
		if(nskip > ncalls) return;
		// calculate restraint energy differences
		// for each residue separately
		for(size_t ir=0;ir< wspace->res.size();ir++){
			double ir_ene = 0;
			for(i = wspace->res[ir].ifirst; i <= wspace->res[ir].ilast; i++) {

				force.diff(wspace->cur.atom[i].p, eqpos[i]);
				dist = force.mag();
				forcemag = 0.5 * k_SI;
				for(int p = 0; p < (Power - 2); p++) {
					forcemag *= dist;
				}
				ene = forcemag;
				ene *= dist * dist;
				ir_ene += ene;
			}
			resene.push_back(ir_ene);
		}

		// now accumulate FEP averages in the usual fashion.

		for(size_t ir=0;ir< wspace->res.size();ir++){
			// First Order terms
			res_a_fep[ir] += exp(-resene[ir]/(Physics::PhysicsConst::kB*300.0));

			for(size_t jr=0;jr< wspace->res.size();jr++){
				// Second Order
				res_ab_fep[ir][jr] += exp(-(resene[ir]+resene[jr])/(Physics::PhysicsConst::kB*300.0));

				for(size_t kr=0;kr< wspace->res.size();kr++){
					// Third Order
					res_abc_fep[ir][jr][kr] += exp(-(resene[ir]+resene[jr]+resene[kr])/
						(Physics::PhysicsConst::kB*300.0));
				}
			}
		}
		ndatapoints++;

		if( (ndatapoints%1000)==0 ){
			printf("t %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f \n",
				-0.592*log(res_a_fep[0]/double(ndatapoints)),
				-0.592*log(res_a_fep[1]/double(ndatapoints)),
				-0.592*log(res_a_fep[2]/double(ndatapoints)),
				-0.592*log(res_a_fep[3]/double(ndatapoints)),
				-0.592*log(res_a_fep[4]/double(ndatapoints)),
				-0.592*log(res_a_fep[5]/double(ndatapoints)),
				-0.592*log(res_a_fep[6]/double(ndatapoints)),
				-0.592*log(res_a_fep[7]/double(ndatapoints)),

				-0.592*log(res_ab_fep[0][1]/double(ndatapoints)),
				-0.592*log(res_ab_fep[0][2]/double(ndatapoints)),
				-0.592*log(res_ab_fep[0][3]/double(ndatapoints)),
				-0.592*log(res_ab_fep[0][4]/double(ndatapoints)),
				-0.592*log(res_ab_fep[0][5]/double(ndatapoints)),
				-0.592*log(res_ab_fep[0][6]/double(ndatapoints)) );
		}
	}

	double FEP_ResConf_Cartesian::calcEnergyAtQ(double newQ){
		return 0.0;
	}

	void FEP_ResConf_Cartesian::printAllData() const{
		double sum_a = 0.0;
		printf("First Order in Matrix: \n");
		for(size_t ir=0;ir< wspace->res.size();ir++){
			printf("one %4d %10.8lf %10.8lf\n",ir,
				res_a_fep[ir]/double(ndatapoints),
				-0.592*log(res_a_fep[ir]/double(ndatapoints)));

			sum_a += -0.592*log(res_a_fep[ir]/double(ndatapoints));
		}
		printf("1st order sum: %10.8f \n",sum_a);

		printf("Second Order in Matrix: \n");
		for(size_t ir=0;ir< wspace->res.size();ir++){
			for(size_t jr=0;jr< wspace->res.size();jr++){
				printf("% 10.8f ",-0.592*log(res_ab_fep[ir][jr]/double(ndatapoints)));
			}
			printf("\n");
		}

		printf("Second Order in Matrix nunbiased:\n");
		for(size_t ir=0;ir< wspace->res.size();ir++){
			for(size_t jr=0;jr< wspace->res.size();jr++){
				printf("%10.8f ",-0.592*log(res_ab_fep[ir][jr]/double(ndatapoints))-
					-0.592*log(res_a_fep[ir]/double(ndatapoints))-
					-0.592*log(res_a_fep[jr]/double(ndatapoints)));
			}
			printf("\n");
		}

		printf("Second Order linear nunbiased:\n");
		for(size_t ir=0;ir< wspace->res.size();ir++){
			for(size_t jr = ir+1;jr< wspace->res.size();jr++){
				printf("two %3d %3d % 10.8f \n", ir, jr,
					-0.592*log(res_ab_fep[ir][jr]/double(ndatapoints))-
					-0.592*log(res_a_fep[ir]/double(ndatapoints))-
					-0.592*log(res_a_fep[jr]/double(ndatapoints)));
			}
		}
		printf("Third Order linear \n");
		for(size_t ir=0;ir< wspace->res.size();ir++){
			for(size_t jr = ir+1;jr< wspace->res.size();jr++){
				for(size_t kr = jr+1;kr< wspace->res.size();kr++){
					printf("threeraw %3d %3d %3d % 10.8f \n", ir, jr, kr,
						-0.592*log(res_abc_fep[ir][jr][kr]/double(ndatapoints)) );
				}
			}
		}
		printf("Third Order linear nunbiased\n");
		for(size_t ir=0;ir< wspace->res.size();ir++){
			for(size_t jr = ir+1;jr< wspace->res.size();jr++){
				for(size_t kr = jr+1;kr< wspace->res.size();kr++){
					printf("three %3d %3d %3d % 10.8f \n", ir, jr, kr,
						-0.592*log(res_abc_fep[ir][jr][kr]/double(ndatapoints))-
						-0.592*log(res_ab_fep[ir][jr]/double(ndatapoints))-
						-0.592*log(res_ab_fep[jr][kr]/double(ndatapoints))-
						-0.592*log(res_ab_fep[ir][kr]/double(ndatapoints))+
						-0.592*log(res_a_fep[ir]/double(ndatapoints))+
						-0.592*log(res_a_fep[jr]/double(ndatapoints))+
						-0.592*log(res_a_fep[kr]/double(ndatapoints)));
				}
			}
		}
	}
}

