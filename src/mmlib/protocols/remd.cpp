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

#include "tools/rdstdout.h"
#include "workspace/workspace.h"
#include "protocols/minimise.h"
#include "protocols/md.h"
#include "forcefields/forcefield.h"

#include "remd.h"

using namespace Physics;
using namespace Maths;
using namespace Monitors;

namespace Protocol{


	REX_Replica::REX_Replica( const Protocol::ProtocolBase &_protocol) : 
		m_Protocol( _protocol.clone() )
	{
			ExchangeMode = Swap;
	}

	void REX_Replica::setInitialSpeeds(){
		int i;
		double sigma;
		double vx, vy, vz;
		
		if( protocol().getWSpace().nAtoms() != m_Structure.nAtoms() ){
			THROW(CodeException,"Error in REX_Replica : Internal SnapShot and WorkSpace have different number of atoms!");
		}

		for(i = 0; i <= m_Structure.nAtoms() ; i++) {
			sigma = sqrt(PhysicsConst::kB * m_Temperature / protocol().getWSpace().atom[i].mass);
			nrand(vx, vy, sigma);
			nrand(vy, vz, sigma);
			m_Structure.atom[i].v.setTo(vx, vy, vz);
		}
	}








	int REX_Local::runcore(){
		int i; // counts over atoms
		int r; // counts over replicas
		int excnt = 0; // counts exchanges
		int Step;
		SnapShot swapm_Structure;
		swapm_Structure = getWSpace().save();
		
		doInternalSafetyCheck();

		double *accfreq = new double[rep.size()];
		int *Tempstat = new int[rep.size() * (Steps+10)]; // 10 for good measure ;)

		for(r = 0; r < rep.size(); r++) {
			accfreq[r] = 0;
			//rep[r].protocol().Steps = ExchangeFreq;// MD runs are each as long as the exchange interval
			//rep[r].protocol().UpdateScr = UpdateScr; // Transfer the local settings to the child processes
			//rep[r].protocol().UpdateTra = UpdateTra; // only the
			//rep[r].protocol().RandVel = false; // MD runs should not randomize the velocities before each run
			rep[r].protocol().OutputLevel = OutputLevel; // If the Parent simulation is silent children are silent too

			if(r == FocusRep) {
				rep[r].protocol().UpdateTra = UpdateTra; // only the focus replic should save trajectory entries
				if(mon.size()>0){
					rep[r].protocol().mon.clear();  // clear any monitors in the template
					rep[r].protocol().mon.add(mon); // the focus rep gets the monitors of the remd
					rep[r].protocol().UpdateMon = UpdateMon;
				}
			} else {
				rep[r].protocol().mon.clear();  // clear any monitors in the template (no monitors at higher temperatures!)
				rep[r].protocol().UpdateTra = -1;
			}

			Tempstat[sqrmat(0, r, rep.size())] = r;
		}

		// Print a header
		if((OutputLevel) &&
			(UpdateScr > 0)){
				printf("remdstart ");
				rep[FocusRep].protocol().infoLineHeader();
				printf("\n");
		}

		// Start the simulation
		for(Step = 0; Step < Steps; Step ++) {
			for(r = 0; r < rep.size(); r++) {
				
				printf("Round: %5d   Replica:  %3d   \n", Step, r);
				
				getWSpace().load(rep[r].m_Structure);
				//rep[r].protocol().m_StepOffset = Step * rep.size();
				//if(Step==0)rep[r].protocol().setup(); 
				
				int result;
				// make sure the temperature is set right
				rep[r].protocol().setTargetTemp( rep[r].getTargetTemp() ); 
				if( Step == 0 )	result = rep[r].protocol().run();
				else            result = rep[r].protocol().runcore();

				if( result < 0 ){
					// if run_core() fails it mean the simulation was unstable and "exploded"
					// in that case attempt to rescue the situation by taking the structure
					// the replica above and copying it with new velocities.
					// Note: this is thermodynamically 'illegal' and can only be tolerated without
					// affecting the results if it occurs sporadically!

					printf("remd.run_core_failure: Temp: %lf \n", rep[r].getTargetTemp());
					int exchangefor = r + 1;
					if(exchangefor >= rep.size())
						exchangefor = r - 1;
					getWSpace().load(rep[exchangefor].m_Structure);

					double factor = sqrt(rep[exchangefor].getTargetTemp() /
						rep[r].getTargetTemp());

					// this replaces the explicit loop below:
					getWSpace().scaleVelocities(1/factor);

					getWSpace().load(rep[r].m_Structure);
					Minimisation mini(*ff);
					mini.Algorithm = Minimisation::ConjugateGradients;
					mini.Steps = 200;
					mini.StepSize = 2E6;
					mini.runcore();
					rep[r].setInitialSpeeds( );
				}

				// if this was the focus replica print an info line
				if(r == FocusRep){
					if((OutputLevel) &&
						(UpdateScr > 0)&&
						 ((excnt % UpdateScr)==0)){
							printf("remd.rep %3d ",excnt);
							rep[r].protocol().infoLine();
							ff->infoLine();
							printf("\n");
					}
				}

				// save structure until next round
				rep[r].m_Structure = getWSpace().save();
			}


			// all replicas have been run. Now attempt exchanges

			// prepare statistics
			for(r = 0; r < rep.size(); r++) {
				Tempstat[sqrmat(excnt + 1, r, rep.size())] = Tempstat[sqrmat(excnt, r, rep.size())];
			}

			// alternatively attempt to swap 0&1, 2&3 , 4&5 ... or 1&2, 3&4, 5&6 ...
			for(r = (excnt % 2); r < (rep.size() - 1); r += 2) {

				for(int u = 0; u < 1; u++) {

					// decide whether to exchange those two structures, by default using the metropolis criterion to
					// REMD delta from the Sugita et al. paper; however this function can be overloaded in a subclass
					// an thus different algorithms can be created.

					if( exchangeCriterion( rep[r] , rep[r + 1] ) ){

						swapm_Structure = rep[r].m_Structure;
						switch ( rep[r].ExchangeMode ){
							case REX_Replica::Swap:
									rep[r].m_Structure = rep[r + 1].m_Structure;
									rep[r + 1].m_Structure = swapm_Structure;
									break;
							case REX_Replica::Downward:
									rep[r].m_Structure = rep[r + 1].m_Structure;
									break;
							case REX_Replica::Upward:
									rep[r + 1].m_Structure = rep[r].m_Structure;
						}

						int rTemp = Tempstat[sqrmat(excnt + 1, r, rep.size())];
						Tempstat[sqrmat(excnt + 1, r, rep.size())] = Tempstat[sqrmat(excnt + 1, r + 1, rep.size())];
						Tempstat[sqrmat(excnt + 1, r + 1, rep.size())] = rTemp;

						// rescale velocities by Temperature factor 
						double factor = sqrt(rep[r + 1].getTargetTemp() / rep[r].getTargetTemp());
						for(i = 0; i < getWSpace().atom.size(); i++) {
							rep[r].m_Structure.atom[i].v.x /= factor;
							rep[r].m_Structure.atom[i].v.y /= factor;
							rep[r].m_Structure.atom[i].v.z /= factor;
							rep[r + 1].m_Structure.atom[i].v.x *= factor;
							rep[r + 1].m_Structure.atom[i].v.y *= factor;
							rep[r + 1].m_Structure.atom[i].v.z *= factor;
						}

						if( OutputLevel >= Verbosity::Loud ){
							printf( "Swapped replicas %d and %d \n", r, r+1 );	
						}

						accfreq[r] += 1.0;
						break;
					}
				}
			}

			excnt++;
		}

		printf("REX_Replica Exchange Statistics: \n");
		for(r = 0; r < rep.size(); r++) {
			accfreq[r] /= (double) excnt;
			printf("rep.stat %3d:\t%4.1lf\t%10.3lf\n",
				r,
				rep[r].getTargetTemp(),
				accfreq[r]);
		}

		printf("REX_Replica exchanges throughout run: \n");
		for(int iex = 0; iex < excnt + 1; iex++) {
			for(r = 0; r < rep.size(); r++)
				printf("%2d\t", Tempstat[sqrmat(iex, r, rep.size())]);
			printf("\n");
		}

		FILE *file = NULL;
		if(rdstdout_file != NULL) file = fopen(rdstdout_file, "at");
		else file = stdout;
		if((file != NULL)&&(file != stdout)) fclose(file);

		delete[]accfreq;
		delete[]Tempstat;

		return Step;
	}

	void REX_Local::getReplicas(SnapShot * m_Structure_ext){
		for(int r = 0; r < rep.size(); r++) {
			m_Structure_ext[r] = rep[r].m_Structure;
		}
	}

	void REX_Local::setAllReplicasTo(SnapShot &m_Structure_ext){
		for(int r = 0; r < rep.size(); r++) {
			rep[r].m_Structure = m_Structure_ext;
		}
	}

	void REX_Local::setReplicas(SnapShot * m_Structure_ext){
		for(size_t r = 0; r < rep.size(); r++) 
		{
			rep[r].m_Structure = m_Structure_ext[r];
		}
	}

	void REX_Local::addReplica(
		const Protocol::ProtocolBase &simTemplate, 
		float Temperature,
		REX_Replica::ExchangeModeType exmode
		)
	{
		REX_Replica newrep(simTemplate);
		newrep.m_Structure = getWSpace().save();
		newrep.protocol().setTargetTemp( Temperature ); 
		// for MD-like protocols, set the velocity vectors accordingly
		newrep.setTargetTemp( Temperature );
		newrep.setInitialSpeeds();
		newrep.ExchangeMode = exmode;
		rep.push_back(newrep);
	}

	void REX_Local::addReplicas(
		const Protocol::ProtocolBase &simTemplate, 
		float firstTemp,
        float factor,
		size_t N,
		REX_Replica::ExchangeModeType exmode
		)
	{
		for(size_t i = 0; i < N; i++)
		{
			addReplica(simTemplate, firstTemp, exmode);
			firstTemp *= factor;
		}
	}


	void REX_Local::addReplicas(
		const Protocol::ProtocolBase &simTemplate, 
		std::vector<float> newtemps,
		REX_Replica::ExchangeModeType exmode
		)
	{
		for(size_t i=0;i<newtemps.size();i++){
			addReplica(simTemplate,newtemps[i],exmode);
		}
	}


	void REX_Local::printCheckPointMIME(){
		for(size_t r = 0; r < rep.size(); r++) 
		{
			rep[r].m_Structure.printMIME(std::string("remd") + int2str(r));
		}
	}

	void REX_Local::readCheckPointMIME(const std::string &filename){
		for(size_t r = 0; r < rep.size(); r++) 
		{
			rep[r].m_Structure.readMIME(filename, std::string("remd") + int2str(r));
		}	
	}

	void REX_Local::doInternalSafetyCheck(){
		if( rep.size() == 0 ){
			throw(ProcedureException("No Replicas were created before running the ReplicaExchange protocol"));
		}

		WorkSpace *wspace = & rep[0].protocol().getWSpace();
		for(size_t r = 0; r < rep.size(); r++) 
		{
			if( ( &rep[0].protocol().getWSpace() ) != wspace ){
				throw(ProcedureException("All Replicas must work on the same workspace !") );
			}
		}	
		

	}

	bool REX_Local::exchangeCriterion( const REX_Replica &lowerrep,
	                                   const REX_Replica &upperrep ) const 
	{
		double delta = (double(lowerrep.m_Structure.epot) - double(upperrep.m_Structure.epot)) *
						(1.0 / (PhysicsConst::kB * upperrep.getTargetTemp()) -
						 1.0 / (PhysicsConst::kB * lowerrep.getTargetTemp()));

		return (delta <= 0) || (exp(-delta) > frand()); 
	}

}

