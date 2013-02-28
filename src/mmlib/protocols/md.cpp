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

// D.  Beeman,  Some Multistep  Methods  for  Use in  Molecular
// Dynamics Calculations, J. Comput. Phys., 20, 130-139 (1976)
//
#include "global.h"

#include "workspace/workspace.h"
#include "protocols/minimise.h"
#include "workspace/space.h"

#include "forcefields/forcefield.h"

#include "protocols/md.h"

using namespace Physics;
using namespace Maths;
using namespace Monitors;

namespace Protocol
{

	MolecularDynamics::MolecularDynamics( Physics::Forcefield & _ff, const PickAtomRange& _def ):
		RangedProtocolBase(_ff,_def)
	{
		settodefault();
	}

	MolecularDynamics::MolecularDynamics( Physics::Forcefield & _ff):
		RangedProtocolBase(_ff)
	{
		settodefault();
	}

	MolecularDynamics::~MolecularDynamics()
	{
	}


	/// sets default parameter values
	void MolecularDynamics::settodefault()
	{
		name = "MolecularDynamics";

		Integrator = MolecularDynamics::Beeman;
		Timestep = 1E-15;
		RandVel = true;
		setTargetTemp(300);

		UpdateRemoveTotalMomentum = 200;
	
		Thermostat = NoThermostat;
		BerendsenTau = 100e-15;
		AndersenRate = 0.1;

		Barostat = NoBarostat;
		TargetPressure = 1;   // Bar
		BerendsenPressureTau = 2000e-15;  // Femtoseconds
		Compress = 0.000046;  // 4.6 E-5 Bar

		LangevinOnHydrogens = true;
		FricCoeff = 10E12;

		Steps = 0;
		m_StepOffset = 0;
		UpdateNList = 10;

		m_CurTemp = 0;
		m_CurPress = 0;

		CentreAfterMove = false;
	};

	void MolecularDynamics::info() const
	{
		printf("---- %s ------------------\n", name.c_str());
		printf("UpdateScr            %d\n", UpdateScr);
		printf("UpdateTra            %d\n", UpdateTra);
		printf("UpdateMon            %d\n", UpdateMon);
		printf("Integrator           ");
		switch (Integrator) {
			case Verlet:					printf("Verlet"); break;
			case VelocityVerlet:	printf("Velocity Verlet"); break;
			case Beeman:					printf("Beeman"); break;
			case Langevin:				printf("Langevin"); break;
			default:							printf("Unknown");
		}
		printf("\n");

		printf("Timestep             %e\n", Timestep );
		printf("Steps                %d\n", Steps );
		printf("RandVel              %d\n", RandVel );
		printf("TargetTemp           ");
		TargetTemp->info();
		printf("\n");
		printf("Thermostat           " );
		switch (Thermostat) {
				case NoThermostat: printf("No Thermostat (constant energy, NVE)"); break;
				case Andersen:       printf("Andersen Thermostat, NVT"); break;
				case Berendsen:      printf("Berendsen Thermostat, NVT"); break;
				default: printf("Unknown");
		}
		printf("\n");
		printf("BerendsenTau         %e\n", BerendsenTau );
		printf("AndersenRate         %lf\n", AndersenRate);
		printf("FricCoeff            %e\n", FricCoeff );
		printf("LangevinOnHydrogens  %s\n", LangevinOnHydrogens ? "true" : "false" );
		printf("Barostat             " );
		switch (Barostat) {
				case NoBarostat: printf("No Barostat"); break;
				case BerendsenBaro:  printf("Berendsen Barostat, NPT"); break;
				default: printf("Unknown");
		}
		printf("\n");
		printf("BerendsenPressureTau %e\n", BerendsenPressureTau );

	}

	void MolecularDynamics::printFinalStatistics(size_t steps)
	{
		printf("Timing: %ld sec (%d hrs %d mins %d secs)",
			(endtime - starttime),
			(endtime - starttime)/3600,
			((endtime - starttime)%3600)/60,
			((endtime - starttime)%3600)%60);
		if(steps>0)
		{
			printf("  %8.5f s/step", double(endtime - starttime)/double(steps));
			printf(" %8.2f hrs/ns", 1E6*double(endtime - starttime)/double(steps)/(3600*Timestep*1E15) );
		}
		if((endtime - starttime)>0){
			printf(" %5.1f ns/day", ( double(steps)*(1E-6*Timestep*1E15) ) * ( 3600 * 24 ) / double(endtime - starttime) );
		}
		printf("\n");
	}

	void MolecularDynamics::setup()
	{
		int start = getStartAtom();
		int end = getEndAtom();
		int natom = getNAtoms();

		m_TargetEkin = 3 * natom * PhysicsConst::kB * TargetTemp->get(0) / 2;

		if(RandVel) { // randomise the velocities if required by user parameter
			printf("Randomizing velocities to boltzmann distribution at T=%.1lf K \n",
				TargetTemp->get(0) );
			setInitialSpeeds(TargetTemp->get(0)); // set speeds to a gaussian dist. with Temperature

			calcTotalLinearVelocity(start, end);
			calcTotalAngularVelocity(start, end);
			calcTotalLinearMomentum(start, end);
			calcTotalAngularMomentum(start, end);

			calcOldPositions();
		}
	}

	int MolecularDynamics::runcore()
	{
		setup();

		if(OutputLevel){
			if((UpdateScr>0)&&(UpdateScr<Steps)){
				getWSpace().boundary().info();
				infoLineHeader();
			}
		}

		m_StatEpot.clear();
		m_StatEtot.clear();
		m_StatTemp.clear();
		m_StatPres.clear();
		m_StatDens.clear();
		m_StatVolu.clear();
		m_TotalEnergyDrift_x.clear();
		m_TotalEnergyDrift_y.clear();

		starttime = (int) time(NULL);
		int stepsTaken = run_core();
		endtime = (int) time(NULL);
		if(OutputLevel) {
 			
			printf("Average potential energy:    %8.2lf (%8.2lf)  kcal/mol\n", m_StatEpot.getAv(),m_StatEpot.getStdDev());
			printf("Average total energy:        %8.2lf (%8.2lf)  kcal/mol\n", m_StatEtot.getAv(),m_StatEtot.getStdDev());
			if(m_TotalEnergyDrift_y.size() > 1){	
				double a,b,R;
				//#for(int jj=0;jj<m_TotalEnergyDrift_y.size();jj++){
			  //#		printf("%f   %f \n",m_TotalEnergyDrift_x[jj],m_TotalEnergyDrift_y[jj]);
				//}
				leastSquaresFit(m_TotalEnergyDrift_x,m_TotalEnergyDrift_y,a,b,R);
				printf("Total energy drift:          %8.4lf           kcal/mol/1000 steps\n",1000.0*b);
			}
			printf("Average temperature:         %8.2lf (%8.2lf)  K\n",        m_StatTemp.getAv(),m_StatTemp.getStdDev());
			printf("Average pressure:            %8.2lf (%8.2lf)  Bar\n",			 m_StatPres.getAv(),m_StatPres.getStdDev());
			printf("Average density:             %8.2lf (%8.2lf)  g/ml\n",	   m_StatDens.getAv(),m_StatDens.getStdDev());
			printf("Average volume:	             %8.2lf (%8.2lf)  A^3\n",			 m_StatVolu.getAv(),m_StatVolu.getStdDev());
			printFinalStatistics(Steps);
		}

		return stepsTaken;
	}

	int MolecularDynamics::run_core()
	{
		int start = getStartAtom();
		int end = getEndAtom();
		int driftfreq = Maths::max((unsigned)1,(unsigned)(Steps/1000));
		m_TotalEnergyDrift_y.clear();
		m_TotalEnergyDrift_x.clear();
		for(Step = 0; Step < Steps; Step++)			// Main principal loop for simulation
		{
			// let particle system know about the current Step nr
			// this also ensures that all the forcefields do a full update at the first step! 
			// (and not use stale information when multistepping is turned on!)
			getWSpace().Step = Step; 
			refreshNeighborList();

			ff->calcForces();
			assertStability();			

			// Do some statistical analysis
			runmonitors();
			m_StatEpot.addPoint(getWSpace().ene.epot * PhysicsConst::J2kcal * PhysicsConst::Na); 
			m_StatTemp.addPoint(m_CurTemp); 
			m_StatPres.addPoint(m_CurPress); 
			m_StatDens.addPoint(getWSpace().getDensity()); 
			m_StatVolu.addPoint(getWSpace().getVolume()); 
 
			applyForces(); // call appropriate Integrator
 //getWSpace().cleanSpace(); // move any stray molecules back into simulation
			if(CentreAfterMove) getWSpace().zeroCentreOfGeometry();
			// add up total energy
			getWSpace().ene.etot += getWSpace().ene.epot;
			getWSpace().ene.etot += getWSpace().ene.ekin;

			m_StatEtot.addPoint(getWSpace().ene.etot * PhysicsConst::J2kcal * PhysicsConst::Na); 
			if(Step > 50)
			if(every(Step,driftfreq)){
				m_TotalEnergyDrift_y.push_back(getWSpace().ene.etot * PhysicsConst::J2kcal * PhysicsConst::Na);
				m_TotalEnergyDrift_x.push_back(Step);
			}
 

			// Analyse (and if nessessary regulate) total translational/rotational momentum
			// this task could maybe be solved more elgantly using a virtual call to the
			// space itself.
			if(every(Step,UpdateRemoveTotalMomentum))
			{
				calcTotalLinearVelocity(start, end);
				calcTotalAngularVelocity(start, end);
				calcTotalLinearMomentum(start, end);
				calcTotalAngularMomentum(start, end);
				if( (Thermostat != MolecularDynamics::NoThermostat) )
				{
					zeroTotalLinearMomentum(start, end);
					zeroTotalAngularMomentum(start, end);
				}
			}

			// Save the coordinates in the trajectory as required
			if(every(Step + m_StepOffset,UpdateTra)) getWSpace().outtra.append();

			// Display any information as required
			if((OutputLevel)&&every(Step + m_StepOffset,UpdateScr)) infoLine();
		}

		return Step;
	}

	void MolecularDynamics::infoLine() const
	{
		// prints a line of current energies
		printf("%10.0lf\t%8.1lf\t%8.1lf\t%8.1lf\t",
			(double) (Step + m_StepOffset) * Timestep * 1E15,
			double(getWSpace().ene.epot) * PhysicsConst::J2kcal * PhysicsConst::Na,
			double(getWSpace().ene.ekin) * PhysicsConst::J2kcal * PhysicsConst::Na,
			double(getWSpace().ene.etot) * PhysicsConst::J2kcal * PhysicsConst::Na);

		printf("%7.1lf",m_CurTemp);
		if((Thermostat != NoThermostat)||
		   (Integrator == Langevin)){
			printf("(%.1lf)",TargetTemp->get(double(Step)/double(Steps)));
		}
		if(getWSpace().getVolume()>=0){
				printf(" %9.2lf",m_CurPress);
				if(Barostat != NoBarostat){
					printf(" (%.2lf) %5.0lf  %e",TargetPressure,getWSpace().getVolume(),getWSpace().ene.InternalVirial);
				}
		}
		mon.printCurData();
		ff->infoLine();
		printf("\n");
	}

	void MolecularDynamics::infoLineHeader() const
	{
		// prints the headers for the above function
		printf("%10s\t%8s\t%8s\t%8s\t", "Time/fs", "Epot", "Ekin", "Etotal");
		
		printf("%7s","T/K");
		if((Thermostat != NoThermostat)||
		   (Integrator == Langevin)){
			printf("(%5s)","tgtT");
		}
		if(getWSpace().getVolume()>=0){
				printf(" %9s","P/Bar");
				if(Barostat != NoBarostat){
					printf(" (%s) %s","tgtP","Vol/A^3");
				}
		}		
		
		mon.printHeader();
		ff->infoLineHeader();
		printf("\n");
	}


	void MolecularDynamics::applyForces(){
		switch (Integrator) {
				case Verlet:
					applyForces_VerletIntegration();
					return;
				case VelocityVerlet:
					applyForces_VelocityVerletIntegration();
					return;
				case Beeman:
					applyForces_BeemanIntegration();
					return;
				case Langevin:
					if(LangevinOnHydrogens){
						applyForces_LangevinIntegration();
					}else{
						applyForces_LangevinIntegration_NonHydrogens();
					}
					return;
		}
	}

	void MolecularDynamics::setInitialSpeeds(){
		setInitialSpeeds(TargetTemp->get(0));
	}

	void MolecularDynamics::setInitialSpeeds(double tgtTemp){
		int i;
		double sigma;
		double vx, vy, vz;
		int natom = getNAtoms();

		if(OutputLevel)
			printf("Setting initital velocities of ensemble.. Temperature %6.1lf\n", tgtTemp);

		int end = getEndAtom();
		for(i = getStartAtom(); i <= end; i++) {
			sigma = sqrt(PhysicsConst::kB * tgtTemp / getWSpace().atom[i].mass);
			nrand(vx, vy, sigma);
			nrand(vy, vz, sigma);
			getWSpace().cur.atom[i].v.setTo(vx, vy, vz);
		}

		calcKineticEnergy();

		setKineticEnergyTo(3 * natom * PhysicsConst::kB * tgtTemp / 2);
		calcKineticEnergy();
	}


	void MolecularDynamics::calcTotalLinearVelocity(int istart, int iend){
		// Set up proxies to workspace to make code more readable.
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array

		int i;
		linvel.setTo(0, 0, 0);
		for(i = istart; i <= iend; i++)
		{
			linvel.add(atom[i].v);
		}
	}

	void MolecularDynamics::calcTotalAngularVelocity(int istart, int iend){
		// set up proxies to workspace to make code more readable.
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array

		dvector p, Li;
		angvel.setTo(0, 0, 0);;
		getWSpace().ene.ekin_angular = 0;

		for(int i = istart; i <= iend; i++) 
		{
			p.setTo(atom[i].p);
			p.mul(PhysicsConst::Angstrom);
			Li.crossProduct(p, atom[i].v);
			angvel.add(Li);
		}

	}

	void MolecularDynamics::calcTotalLinearMomentum(int istart, int iend){
		// set up proxies to workspace to make code more readable.
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array

		dvector p;
		linmom.setTo(0, 0, 0);
		double totmass = 0;

		for(int i = istart; i <= iend; i++) 
		{
			p.setTo(atom[i].v);
			p.mul(getWSpace().atom[i].mass);
			totmass += getWSpace().atom[i].mass;
			linmom.add(p);
		}
		getWSpace().ene.ekin_linear = 0.5 * sqr(linmom.mag()) / totmass;
	}

	void MolecularDynamics::calcTotalAngularMomentum(int istart, int iend)
	{
		// set up proxies to workspace to make code more readable.
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array

		int i;
		dvector p, Li;
		double myL = 0, myLi;
		angmom.setTo(0, 0, 0);
		double Inertia = 0;
		double totmass = 0;

		dvector com; // centre of mass

		//first calculate centre of mass
		com.setTo(0, 0, 0);
		for(i = istart; i <= iend; i++) 
		{
			p.setTo(atom[i].p);
			p.mul(getWSpace().atom[i].mass);
			totmass += getWSpace().atom[i].mass;
			com.add(p);
		}
		com.div(totmass);

		for(i = istart; i <= iend; i++) 
		{
			p.setTo(atom[i].p);
			p.sub(com);
			p.mul(PhysicsConst::Angstrom);
			Li.crossProduct(p, atom[i].v, getWSpace().atom[i].mass);

			Inertia += getWSpace().atom[i].mass * dotProduct(p, p);

			myLi = (getWSpace().atom[i].mass * atom[i].v.mag() * p.mag());
			myL += myLi;

			//printf("myLi %e Li.x/Li.y/Li.z %e %e %e |Li| %e \n", myLi, Li.x, Li.y, Li.z, Li.mag());
			angmom.add(Li);
		}
		//printf("scalar L %e\n",myL);
		//printf("dvector L %e\n",angmom.mag());

		getWSpace().ene.ekin_angular = 0.5 * sqr((double) angmom.mag()) / Inertia;
	}


	void MolecularDynamics::zeroTotalLinearMomentum(int istart, int iend)
	{
		// set up proxies to workspace to make code more readable.
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array

		int i;
		dvector ptot;

		for(i = istart; i <= iend; i++) 
		{
			ptot.setTo(linmom);
			ptot.div((iend - istart) * getWSpace().atom[i].mass);
			atom[i].v.sub(ptot);
		}

	}

	void MolecularDynamics::zeroTotalAngularMomentum(int istart, int iend){
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array
		int i;
		dvector p, ptot;

		dvector com;
		com = getWSpace().getCentreOfMass();

		for(i = istart; i <= iend; i++) 
		{
			p.setTo(atom[i].p);
			p.sub(com);
			p.mul(PhysicsConst::Angstrom);

			ptot.crossProduct(angmom, p);
			ptot.unify();
			ptot.mul(angmom.mag() / p.mag());
			ptot.div((iend - istart) * getWSpace().atom[i].mass);

			//printf(" %e \n",ptot.mag());
			atom[i].v.sub(ptot);
		}

	}


	void MolecularDynamics::calcKineticEnergy()
	{
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array
		int natom = getNAtoms();
		int i;
		getWSpace().ene.ekin = 0;

		int end = getEndAtom();
		for(i = getStartAtom(); i <= end; i++)
		{
			getWSpace().ene.ekin += (0.5 * getWSpace().atom[i].mass * sqr(atom[i].v.mag()));
		}

		// Calculate Temperature from K = 2NkBT/2
		m_CurTemp = 2.0 * getWSpace().ene.ekin / (3 * natom * PhysicsConst::kB);

		// Get volume and convert into units of meters (instead of Angstroms)
		double V = getWSpace().getVolume() * 1E-30;
		// V will be < 0 if it's an infinite space and thus no pressure can be calculated
		if(V > 0.0){
			// Calculate Pressure from P = (2.0K + W)/(3V)
			// Where W is the internal virial which is basically the derivative of the internal 
			// energy wrt to the volume, i.e. dH/dV
			//
			//  dH/dV = W/(3.0*V)
			//  For pairwise interactions W is: (see virial calculation in the forcefield components)
			//  W = Sum[ Dist * |Force|  ] 
			//
			// The first term  (2.0K/3V) is basically the ideal gas contribution PV = nRT.

			m_CurPress = (2.0*getWSpace().ene.ekin - getWSpace().ene.InternalVirial)/(3.0 * V) / PhysicsConst::Bar2Pa;


		}

	}


	// adjusts the speeds of all the particles such that the total
	// kinetic energy equals the desired value; also sets the velocities so that the
	// total momentum = 0;

	void MolecularDynamics::setKineticEnergyTo(double Edes)
	{
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array
		double chi;
		int i;

		chi = sqrt(Edes / double(getWSpace().ene.ekin));

		int end = getEndAtom();
		for(i = getStartAtom(); i <= end; i++)
		{
			atom[i].v.mul(chi);
		}
	}


	// Applies Thermostat as defined in config file
	void MolecularDynamics::applyThermostat()
	{
		switch (Thermostat) {
			case NoThermostat:
				return;
			case Andersen:
				applyAndersonThermostat(TargetTemp->get(double(Step)/double(Steps)));
				return;
			case Berendsen:
				applyBerendsenThermostat(TargetTemp->get(double(Step)/double(Steps)));
				return;
		}
	}

	// adjusts the speeds of all the particles such that the total
	// kinetic energy equals the desired value more;

	void MolecularDynamics::applyBerendsenThermostat(double TargetTargetTemp)
	{
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array
		int natom = getNAtoms();
		double chi; // adjustment factor
		double Edes; // desired Ekinetic
		int i;

		Edes = 3.0 * natom * PhysicsConst::kB * TargetTargetTemp / 2.0;
		chi = sqrt(1.0 + Timestep * (Edes / double(getWSpace().ene.ekin) - 1.0) / BerendsenTau);

		int end = getEndAtom();
		for(i = getStartAtom(); i <= end; i++)
		{
			atom[i].v.mul(chi);
		}
	}


	void MolecularDynamics::applyAndersonThermostat(double TargetTargetTemp)
	{
		int i;
		double sigma;
		double vx, vy, vz;
		double rate;
		int intrate;
		dvector vt, t2a;
		int natom = getNAtoms();

		// change as many atoms such that natom atoms are changed after Temperature_tau has passed
		// note that this is very naive and may well have to be replaced

		rate = AndersenRate * Timestep * 1E15 / pow(double(natom), double( 2.0 / 3.0 ));

		intrate = (int) (10000.0 * rate);

		int end = getEndAtom();
		for(i = getStartAtom(); i <= end; i++)
		{
			// randomly decide if this atom is to be changed.
			if((rand() % 10000) > intrate)
				continue;

			//sample a new velocity from a Maxwell-Boltzmann distribution
			sigma = sqrt(PhysicsConst::kB * TargetTargetTemp / getWSpace().atom[i].mass);
			nrand(vx, vy, sigma);
			nrand(vy, vz, sigma);

			// avoid exremely fast particles (tends to crash the simulation)
			if(fabs(vx) > (sigma * 4))
				continue;
			if(fabs(vy) > (sigma * 4))
				continue;
			if(fabs(vz) > (sigma * 4))
				continue;

			getWSpace().cur.atom[i].v.setTo(vx, vy, vz);

			//printf("%8.3lf\n", sqrt(vx*vx + vy*vy + vz*vz));
			//printf("%8.3lf\n%8.3lf\n%8.3lf\n", vx,vy,vz);
		}
	}



	// Applies Thermostat as defined in config file
	void MolecularDynamics::applyBarostat()
	{
	 
		switch (Barostat) {
			case NoBarostat:
				return;
			case BerendsenBaro:
				applyBerendsenBarostat();
				return;
		}
		 
	}

	
	void MolecularDynamics::applyBerendsenBarostat()
	{
		double factor;  //scaling factor
		factor = pow(1.0 + (Timestep*Compress*
			                 (m_CurPress - TargetPressure)
											 /(BerendsenPressureTau)),1.0/3.0);
		getWSpace().scaleSystem(factor);
	}


	// this is to prevent hot spots
	void MolecularDynamics::equaliseVelocities()
	{
		int natom = getNAtoms();
		int *index = new int[natom];
		double *iv = new double[natom];
		int arraysize = 0;
		int i;

		int end = getEndAtom();
		for(i = getStartAtom(); i <= end; i++)
		{
			index[arraysize] = i;
			iv[arraysize] = getWSpace().cur.atom[i].v.mag();
			arraysize++;
		}

		qcksort(&iv[0], &index[0], arraysize);

		//for(i = 0; i < (arraysize / 2); i++) {
		//	// factor = (sqr(iv[index[i]])+sqr([index[arraysize-i-1]]))
		//}

		delete[]iv;
		delete[]index;
	}

	void MolecularDynamics::calcOldPositions()
	{
		// Set up proxies to workspace to make code more readable.
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array
		SnapShotAtom *oldatom = getWSpace().old.atom; // atom coordinate array

		int i;
		dvector vt, t2a;

		int end = getEndAtom();
		for(i = getStartAtom(); i <= end; i++)
		{
			t2a.setTo(atom[i].f);
			t2a.setTo(atom[i].f);
			t2a.mul(0.5 * sqr(Timestep) / PhysicsConst::Angstrom / getWSpace().atom[i].mass);

			vt.setTo(atom[i].v);
			vt.mul(Timestep / PhysicsConst::Angstrom);

			oldatom[i].p.setTo(atom[i].p);
			oldatom[i].p.sub(vt);
			oldatom[i].p.add(t2a);
		}
	}

	void MolecularDynamics::reCalculateNewPositions()
	{
		// Set up proxies to workspace to make code more readable.
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array
		SnapShotAtom *oldatom = getWSpace().old.atom; // atom coordinate array

		int i;
		dvector vt;

		int end = getEndAtom();
		for(i = getStartAtom(); i <= end; i++)
		{
			vt.setTo(atom[i].v);
			vt.mul(Timestep / PhysicsConst::Angstrom);

			atom[i].p.setTo(oldatom[i].p);
			atom[i].p.add(vt);
		}
	}


	void MolecularDynamics::applyForces_VerletIntegration()
	{
		// Set up proxies to workspace to make code more readable.
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array
		SnapShotAtom *oldatom = getWSpace().old.atom; // atom coordinate array

		int i;
		dvector np,t2a;
		dvector dv;
		dvector dr, drverlet;
		dvector oldv;
		
		int end = getEndAtom();
		for(i = getStartAtom(); i <= end; i++)
		{
			t2a.setTo(atom[i].f);
			t2a.mul((sqr(Timestep) / PhysicsConst::Angstrom / getWSpace().atom[i].mass));

			np.setTo(atom[i].p);
			np.mul(2);
			np.sub(oldatom[i].p);
			np.add(t2a); // add all velvet terms to get new postiions

			atom[i].v.setTo(np);
			atom[i].v.sub(oldatom[i].p);
      getWSpace().boundary().getClosestImage(atom[i].v);
			atom[i].v.div(2 * Timestep / PhysicsConst::Angstrom); // calculate velocities ( before loosing r(t-dt) )

			oldatom[i].p.setTo(atom[i].p);
			atom[i].p.setTo(np);
		}

		//remember old forces
		for(i = getStartAtom(); i <= end; i++)
		{
			oldatom[i].f.setTo(atom[i].f);
		}

		calcKineticEnergy();
		applyThermostat();   //adjust velocities according to Thermostat
		applyBarostat();     //scale the system according to Barostat

	}


	void MolecularDynamics::applyForces_VelocityVerletIntegration(){
		// Set up proxies to workspace to make code more readable.
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array

		int i;
		dvector dr, dv;
		double t = Timestep;
		double tinvmass;
		getWSpace().ene.ekin = 0;

		int end = getEndAtom();
		for(i = getStartAtom(); i <= end; i++)
		{
			tinvmass = 0.5 * t / getWSpace().atom[i].mass; // t div m

			// finish change in velocity using a(n+1)
			dv.x = atom[i].f.x * tinvmass;
			dv.y = atom[i].f.y * tinvmass;
			dv.z = atom[i].f.z * tinvmass;
			atom[i].v.add(dv);
		}

		// The velocities are now "real", i.e. finished and thus now is the time
		// to calculate the kinetic energy and apply any barostat
		// note no other thermostats are called here because 
		// langevin dynamics regultes the temperature itself - thermostat
		// is built-in so to speak
		calcKineticEnergy();
		applyThermostat();   //adjust velocities according to Thermostat
		applyBarostat();     //scale the system according to Barostat

		for(i = getStartAtom(); i <= end; i++)
		{
			tinvmass = 0.5 * t / getWSpace().atom[i].mass; // t div m

			// finish change in velocity using a(n+1)
			dv.x = atom[i].f.x * tinvmass;
			dv.y = atom[i].f.y * tinvmass;
			dv.z = atom[i].f.z * tinvmass;
			atom[i].v.add(dv);

			// change in position
			dr.x = t * (atom[i].v.x + dv.x);
			dr.y = t * (atom[i].v.y + dv.y);
			dr.z = t * (atom[i].v.z + dv.z);
			dr.mul(1 / PhysicsConst::Angstrom);

			atom[i].p.add(dr);
			atom[i].v.add(dv);
		}
	}

	// MolecularDynamics::Beeman's algorithm
	void MolecularDynamics::applyForces_BeemanIntegration(){
		// Set up proxies to workspace to make code more readable.
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array
		SnapShotAtom *oldatom = getWSpace().old.atom; // atom coordinate array
		int start = getStartAtom();
		int end = getEndAtom();

		int i;
		double t;
		double t2;
		double tinvmass;
		dvector dr, dv;
		dvector t2a;

		getWSpace().ene.ekin = 0;

		t = Timestep;
		t2 = sqr(t);

		if(Step > 0) {  // dont do this on first Step
			for(i = start; i <= end; i++)
			{
				tinvmass = t / (6.0 * getWSpace().atom[i].mass); // t div 6m

				// finish change in velocity using a(n+1)
				dv.x = ((2.0) * atom[i].f.x) * tinvmass;
				dv.y = ((2.0) * atom[i].f.y) * tinvmass;
				dv.z = ((2.0) * atom[i].f.z) * tinvmass;
				atom[i].v.add(dv);
			}
		}

		calcKineticEnergy(); //calculate instanteneous temperature & pressure
		applyThermostat();   //adjust velocities according to Thermostat
		applyBarostat();     //scale the system according to Barostat

		if(Step < (Steps - 1)) {
			// predict new positions
			for(i = start; i <= end; i++)
			{
				tinvmass = t / (6.0 * getWSpace().atom[i].mass); // t div 6m

				// change in position
				dr.x = t * (atom[i].v.x + (4.0 * atom[i].f.x - oldatom[i].f.x) * tinvmass);
				dr.y = t * (atom[i].v.y + (4.0 * atom[i].f.y - oldatom[i].f.y) * tinvmass);
				dr.z = t * (atom[i].v.z + (4.0 * atom[i].f.z - oldatom[i].f.z) * tinvmass);

				dr.mul(1 / PhysicsConst::Angstrom);

				dv.x = (5.0 * atom[i].f.x - oldatom[i].f.x) * tinvmass;
				dv.y = (5.0 * atom[i].f.y - oldatom[i].f.y) * tinvmass;
				dv.z = (5.0 * atom[i].f.z - oldatom[i].f.z) * tinvmass;

				//update positions and speeds.
				atom[i].p.add(dr);
				atom[i].v.add(dv);

				// save current force as next Step's old force
				oldatom[i].f.setTo(atom[i].f);
			}
		}
	}



	void MolecularDynamics::applyForces_LangevinIntegration(){
		// Set up proxies to workspace to make code more readable.
		SnapShotAtom *atom = getWSpace().cur.atom; // atom coordinate array
		int start = getStartAtom();
		int end = getEndAtom();

		int i;
		dvector dr, dv;
		double t = Timestep;
		double tinvmass;
		getWSpace().ene.ekin = 0;

		double gamma = FricCoeff;
		double gammat = gamma * t;
		double gt = gammat; // these are the consecutive powers for the series expnsions
		double gt2 = gt * gt;
		double gt3 = gt * gt2;
		double gt4 = gt2 * gt2;
		double gt5 = gt2 * gt3;
		double gt6 = gt3 * gt3;
		double gt7 = gt3 * gt4;
		double gt8 = gt4 * gt4;
		double gt9 = gt4 * gt5;

		double c0; // langevin coefficients
		double c1;
		double c2;

		double T = TargetTemp->get(double(Step)/double(Steps)); // random force calculation
		double egammat;
		double varv;
		double varr;
		double crv;

		double kTdivm;
		double sigmav;
		double sigmar;
		double n1x, n2x;
		double n1y, n2y;
		double n1z, n2z;
		dvector deltav;
		dvector deltar;

		if(gammat < 0.00001) {
			// when the friction coeffcient is really small langevin dynamics
			// approximates simple velocity verlet integration so just call that
			// instead (its faster)
			applyForces_VelocityVerletIntegration();
			return;
		} else if(gammat > 0.05) {
			// when the friction coeffcient is relatively large, run langevin dynamics
			// the parameters must be calculated using exponential functions
			c0 = exp(-gammat); // langevin coefficients
			c1 = (1 - c0) / gammat;
			c2 = (1 - c1) / gammat;

			egammat = exp(-gammat);
			varv = (1.0 - sqr(egammat));
			varr = (2.0 * gammat - 3.0 + (4.0 - egammat) * egammat);
			crv = sqr(1.0 - egammat) / (sqrt(varv * varr));
		} else {
			// when the friction coeffcient is in the midrange we can approximate the paramters
			// by taylor expansion which is a little faster
			c0 =
				1 - gt + 0.5 * gt2 - gt3 / 6.0 + gt4 / 24.0 - gt5 / 120.0 + gt6 / 720.0 - gt7 / 5040.0 + gt8 / 40320.0 -
				gt9 / 362880.0;
			c1 =
				1 - (0.5 * gt - gt2 / 6.0 + gt3 / 24.0 - gt4 / 120.0 + gt5 / 720.0 - gt6 / 5040.0 + gt7 / 40320.0 -
				gt8 / 362880.0);
			c2 =
				0.5 - gt / 6.0 + gt2 / 24.0 - gt3 / 120.0 + gt4 / 720.0 - gt5 / 5040.0 + gt6 / 40320.0 - gt7 / 362880.0;

			varr = 2.0 * gt3 / 3.0 - gt4 / 2.0 + 7.0 * gt5 / 30.0 - gt6 / 12.0
				+ 31.0 * gt7 / 1260.0 - gt8 / 160.0 + 127.0 * gt9 / 90720.0;
			varv = 2.0 * gt - 2.0 * gt2 + 4.0 * gt3 / 3.0 - 2.0 * gt4 / 3.0 + 4.0 * gt5 / 15.0
				- 4.0 * gt6 / 45.0 + 8.0 * gt7 / 315.0 - 2.0 * gt8 / 315.0 + 4.0 * gt9 / 2835.0;
			crv = sqrt(3.0) * (0.5 - 3.0 * gt / 16.0 - 17.0 * gt2 / 1280.0 + 17.0 * gt3 / 6144.0
				+ 40967.0 * gt4 / 34406400.0 - 57203.0 * gt5 / 275251200.0 -
				1429487.0 * gt6 / 13212057600.0);
		}

		if(Step > 0) {
			for(i = start; i <= end; i++)
			{
				tinvmass = t / getWSpace().atom[i].mass; // t div m

				// finish change in velocity using a(n+1)
				dv.x = (1.0 - c0 / c1) * atom[i].f.x * tinvmass / gammat;
				dv.y = (1.0 - c0 / c1) * atom[i].f.y * tinvmass / gammat;
				dv.z = (1.0 - c0 / c1) * atom[i].f.z * tinvmass / gammat;
				atom[i].v.add(dv);
			}
		}

		// The velocities are now "real", i.e. finished and thus now is the time
		// to calculate the kinetic energy and apply any barostat
		// note no other thermostats are called here because 
		// langevin dynamics regultes the temperature itself - thermostat
		// is built-in so to speak
		calcKineticEnergy();
		applyBarostat();     //scale the system according to Barostat

		if(Step < (Steps - 1)) {
			for(i = start; i <= end; i++)
			{
				tinvmass = t / getWSpace().atom[i].mass; // t div m

				kTdivm = PhysicsConst::kB * T / getWSpace().atom[i].mass;
				sigmav = sqrt(kTdivm * varv);
				sigmar = sqrt(kTdivm * varr) / gamma;

				nrand(n1x, n2y, (double) 1);
				nrand(n2x, n1z, (double) 1);
				nrand(n1y, n2z, (double) 1);

				deltav.setTo(sigmav * (crv * n1x + sqrt(1.0 - sqr(crv)) * n2x),
					sigmav * (crv * n1y + sqrt(1.0 - sqr(crv)) * n2y),
					sigmav * (crv * n1z + sqrt(1.0 - sqr(crv)) * n2z));
				deltar.setTo(sigmar * n1x, sigmar * n1y, sigmar * n1z);

				// change in position
				dr.x = t * (c1 * atom[i].v.x + c2 * atom[i].f.x * tinvmass) + deltar.x;
				dr.y = t * (c1 * atom[i].v.y + c2 * atom[i].f.y * tinvmass) + deltar.y;
				dr.z = t * (c1 * atom[i].v.z + c2 * atom[i].f.z * tinvmass) + deltar.z;
				dr.mul(1 / PhysicsConst::Angstrom);

				// finish change in velocity using a(n+1)
				dv.x = (c0 * c2) * atom[i].f.x * tinvmass / c1 + deltav.x;
				dv.y = (c0 * c2) * atom[i].f.y * tinvmass / c1 + deltav.y;
				dv.z = (c0 * c2) * atom[i].f.z * tinvmass / c1 + deltav.z;

				atom[i].p.add(dr);
				atom[i].v.mul(c0);
				atom[i].v.add(dv);
			}
		}
	}

	
	/// A specialized version of the normal Langevin integration that 
	/// restricts  stochastic behaviour to 
	/// non hydrogren atoms only - the remaining atoms are integrated using
	/// velocity verlet integration (i.e. strict newtonian mechanics)
	void MolecularDynamics::applyForces_LangevinIntegration_NonHydrogens(){
		// Set up proxies to workspace to make code more readable.
		SnapShotAtom *atom = getWSpace().cur.atom;   // atom coordinate array
		ParticleStore& atomparam = getWSpace().atom; // atom property array

		int start = getStartAtom();
		int end = getEndAtom();

		int i;
		dvector dr, dv;
		double t = Timestep;
		double tinvmass;
		getWSpace().ene.ekin = 0;

		double gamma = FricCoeff;
		double gammat = gamma * t;
		double gt = gammat; // these are the consecutive powers for the series expnsions
		double gt2 = gt * gt;
		double gt3 = gt * gt2;
		double gt4 = gt2 * gt2;
		double gt5 = gt2 * gt3;
		double gt6 = gt3 * gt3;
		double gt7 = gt3 * gt4;
		double gt8 = gt4 * gt4;
		double gt9 = gt4 * gt5;

		double c0; // langevin coefficients
		double c1;
		double c2;

		double T = TargetTemp->get(double(Step)/double(Steps)); // random force calculation
		double egammat;
		double varv;
		double varr;
		double crv;

		double kTdivm;
		double sigmav;
		double sigmar;
		double n1x, n2x;
		double n1y, n2y;
		double n1z, n2z;
		dvector deltav;
		dvector deltar;

		if(gammat < 0.00001) {
			// when the friction coeffcient is really small langevin dynamics
			// approximates simple velocity verlet integration so just call that
			// instead (its faster)
			applyForces_VelocityVerletIntegration();
			return;
		} else if(gammat > 0.05) {
			// when the friction coeffcient is relatively large, run langevin dynamics
			// the parameters must be calculated using exponential functions
			c0 = exp(-gammat); // langevin coefficients
			c1 = (1 - c0) / gammat;
			c2 = (1 - c1) / gammat;

			egammat = exp(-gammat);
			varv = (1.0 - sqr(egammat));
			varr = (2.0 * gammat - 3.0 + (4.0 - egammat) * egammat);
			crv = sqr(1.0 - egammat) / (sqrt(varv * varr));
		} else {
			// when the friction coeffcient is in the midrange we can approximate the paramters
			// by taylor expansion which is a little faster
			c0 =
				1 - gt + 0.5 * gt2 - gt3 / 6.0 + gt4 / 24.0 - gt5 / 120.0 + gt6 / 720.0 - gt7 / 5040.0 + gt8 / 40320.0 -
				gt9 / 362880.0;
			c1 =
				1 - (0.5 * gt - gt2 / 6.0 + gt3 / 24.0 - gt4 / 120.0 + gt5 / 720.0 - gt6 / 5040.0 + gt7 / 40320.0 -
				gt8 / 362880.0);
			c2 =
				0.5 - gt / 6.0 + gt2 / 24.0 - gt3 / 120.0 + gt4 / 720.0 - gt5 / 5040.0 + gt6 / 40320.0 - gt7 / 362880.0;

			varr = 2.0 * gt3 / 3.0 - gt4 / 2.0 + 7.0 * gt5 / 30.0 - gt6 / 12.0
				+ 31.0 * gt7 / 1260.0 - gt8 / 160.0 + 127.0 * gt9 / 90720.0;
			varv = 2.0 * gt - 2.0 * gt2 + 4.0 * gt3 / 3.0 - 2.0 * gt4 / 3.0 + 4.0 * gt5 / 15.0
				- 4.0 * gt6 / 45.0 + 8.0 * gt7 / 315.0 - 2.0 * gt8 / 315.0 + 4.0 * gt9 / 2835.0;
			crv = sqrt(3.0) * (0.5 - 3.0 * gt / 16.0 - 17.0 * gt2 / 1280.0 + 17.0 * gt3 / 6144.0
				+ 40967.0 * gt4 / 34406400.0 - 57203.0 * gt5 / 275251200.0 -
				1429487.0 * gt6 / 13212057600.0);
		}

		if(Step > 0) {
			for(i = start; i <= end; i++)
			{
				if(atomparam[i].isHydrogen()){
					tinvmass = 0.5 * t / getWSpace().atom[i].mass; // t div m

					// finish change in velocity using a(n+1)
					dv.x = atom[i].f.x * tinvmass;
					dv.y = atom[i].f.y * tinvmass;
					dv.z = atom[i].f.z * tinvmass;
					atom[i].v.add(dv);
				}
				else
				{
					tinvmass = t / getWSpace().atom[i].mass; // t div m

					// finish change in velocity using a(n+1)
					dv.x = (1.0 - c0 / c1) * atom[i].f.x * tinvmass / gammat;
					dv.y = (1.0 - c0 / c1) * atom[i].f.y * tinvmass / gammat;
					dv.z = (1.0 - c0 / c1) * atom[i].f.z * tinvmass / gammat;
					atom[i].v.add(dv);
				}
			}
		}

		// The velocities are now "real", i.e. finished and thus now is the time
		// to calculate the kinetic energy and apply any barostat
		// note no other thermostats are called here because 
		// langevin dynamics regultes the temperature itself - thermostat
		// is built-in so to speak
		calcKineticEnergy();
		applyBarostat();     //scale the system according to Barostat

		if(Step < (Steps - 1)) {
			for(i = start; i <= end; i++)
			{
				if(atomparam[i].isHydrogen()){
					tinvmass = 0.5 * t / getWSpace().atom[i].mass; // t div m

					// half if change in velocity using a(n) (the rest is done in the next step)
					dv.x = atom[i].f.x * tinvmass;
					dv.y = atom[i].f.y * tinvmass;
					dv.z = atom[i].f.z * tinvmass;

					// change in position
					dr.x = t * (atom[i].v.x + dv.x);
					dr.y = t * (atom[i].v.y + dv.y);
					dr.z = t * (atom[i].v.z + dv.z);
					dr.mul(1 / PhysicsConst::Angstrom);

					atom[i].p.add(dr);
					atom[i].v.add(dv);
				}
				else
				{
					tinvmass = t / getWSpace().atom[i].mass; // t div m

					kTdivm = PhysicsConst::kB * T / getWSpace().atom[i].mass;
					sigmav = sqrt(kTdivm * varv);
					sigmar = sqrt(kTdivm * varr) / gamma;

					nrand(n1x, n2y, (double) 1);
					nrand(n2x, n1z, (double) 1);
					nrand(n1y, n2z, (double) 1);

					deltav.setTo(sigmav * (crv * n1x + sqrt(1.0 - sqr(crv)) * n2x),
						sigmav * (crv * n1y + sqrt(1.0 - sqr(crv)) * n2y),
						sigmav * (crv * n1z + sqrt(1.0 - sqr(crv)) * n2z));
					deltar.setTo(sigmar * n1x, sigmar * n1y, sigmar * n1z);

					// change in position
					dr.x = t * (c1 * atom[i].v.x + c2 * atom[i].f.x * tinvmass) + deltar.x;
					dr.y = t * (c1 * atom[i].v.y + c2 * atom[i].f.y * tinvmass) + deltar.y;
					dr.z = t * (c1 * atom[i].v.z + c2 * atom[i].f.z * tinvmass) + deltar.z;
					dr.mul(1 / PhysicsConst::Angstrom);

					// finish change in velocity using a(n+1)
					dv.x = (c0 * c2) * atom[i].f.x * tinvmass / c1 + deltav.x;
					dv.y = (c0 * c2) * atom[i].f.y * tinvmass / c1 + deltav.y;
					dv.z = (c0 * c2) * atom[i].f.z * tinvmass / c1 + deltav.z;

					atom[i].p.add(dr);
					atom[i].v.mul(c0);
					atom[i].v.add(dv);
				}
			}
		}

	}



} // namespace 'Protocol'

