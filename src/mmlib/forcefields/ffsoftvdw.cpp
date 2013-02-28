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

#include "workspace/neighbourlist.h"
#include "workspace/bondorder.h"
#include "forcefields/ffsoftvdw.h"

//namespace includes
using namespace Maths;

namespace Physics
{
	void FF_SoftVDW::setup()
	{
		// Proxies
		WorkSpace& wspace = getWSpace();
		Active = true;
		wspace.nlist().requestCutoff(VdwCutoff);
		if(AutoScaling14)
		{
			Vdw14Scaling = wspace.ffps().Vdw14Scaling;
		}
	}

	void FF_SoftVDW::info() const // prints a little block of parameter information
	{
		printf(" --FF_SoftVDW -------------- \n");
		if( m_Hard )
		{
		}
		else
		{
			printf("SoftVDW Hardness Factor:   %6.3lf\n", Hardness);			
		}
		printf("\n");
	}

	void FF_SoftVDW::infoLineHeader() const // prints the headers for the above function
	{
		printf("\t%6s", "ESoftVDW");
	}

	void FF_SoftVDW::infoLine() const // prints a line of current energies
	{
		// Hamiltonian does not currently contain this component ...
		Hamiltonian *printene = &getWSpace().ene;
		printf("%8.1lf", (double)(printene->epot_vdw * PhysicsConst::J2kcal * PhysicsConst::Na ) );
	}

	void FF_SoftVDW::calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level)
	{
		calcEnergies();
	}

	void FF_SoftVDW::calcEnergies()
	{
		calcForces();
	}

	void FF_SoftVDW::calcForces()
	{
		if( m_Hard )
		{
			calcForces_NormalVDW();
		}
		else
		{
			calcForces_Soft();
		}
	}

	void FF_SoftVDW::calcForces_Soft()
	{
		if(!Active) return;

		// Place this in the correct unit system
		const double localHardness =  Hardness * (PhysicsConst::kcal2J / PhysicsConst::Na);
		const double fourInvAnglocalHardness = 4.0 * PhysicsConst::invAngstrom * localHardness;
		const double sqrVdwCutoff = VdwCutoff * VdwCutoff;

		// set up proxies to workspace to make code more readable.
		WorkSpace& wspace = getWSpace();
		size_t natom = wspace.atom.size();             // number of atoms in workspace
		const ParticleStore& atomparam = wspace.atom;  // atom parameter array
		SnapShotAtom *atom = wspace.cur.atom;						// atom coordinate array
		const NeighbourData *fnbor = wspace.nlist().getData(); // neighborlist

		int i, j, nj; // i&j are atom indices, nj count through neighbor list
		Tvector < double > fv; // force Tvector<double>
		double soft_force = 0.0; // individual force magnitudes
		double soft_potential = 0.0; // individual potential energy contributions

		double vdw_potential = 0; // individual potential energy contributions
		double vdw_force = 0; // individual potential energy contributions
		double atomi_radius;		

		// reset the potential energy
		epot = 0;
		wspace.ene.epot_vdw = 0;
		wspace.ene.epot_vdw_att = 0;
		wspace.ene.epot_vdw_rep = 0;

		for(i = 0; i < natom; i++)
		{ 
			// loop over all particles
			atomi_radius = atomparam[i].radius;

			for(nj = 0; nj < fnbor[i].n; nj++) // and all their neighbours
			{
				j = NList32Bit_Index(fnbor[i].i[nj]);   // get atom index index
				if(j > i) break;                        // no shadow neighbors

				if( NList32Bit_BondOrder(fnbor[i].i[nj]) <= BondOrder_1_4_Pair) continue; // only do atoms with bond order > 1-4

				double sqrdistij = wspace.atomxyz(i).sqrdist( wspace.atomxyz(j) ); // distance in Angstrom, *NOT* SI units !!
				double VDWSum = atomi_radius + atomparam[j].radius;
				VDWSum -= OlapOffset;
				double sqrVDWSum = VDWSum * VDWSum;

				if(sqrdistij < sqrVDWSum && sqrdistij < sqrVdwCutoff) // then we want to evaluate the energy contribution
				{
					double distij = sqrt( sqrdistij );
					double A = 1.0/VDWSum;
					double A2X = A * A * distij;
					double A2X2 = A2X * distij;
					vdw_potential = (1.0 - A2X2);
					vdw_potential *= vdw_potential; // square it
					vdw_potential *= localHardness;
					vdw_force = fourInvAnglocalHardness * ( (A2X2 * A2X) - (A2X) );

					// add up potentials -----------------------------
					wspace.ene.epot += vdw_potential;
					wspace.ene.epot_vdw += vdw_potential;

					// Now apply the force.
					fv.setTo(atom[i].p);
					fv.sub(atom[j].p); // fv is now the vector along which all forces will act. fv points *from* atom i *to* atom nj
					fv.mul(vdw_force / distij);
					atom[i].f.sub(fv);
					atom[j].f.add(fv);
				}
			}
		}

		epot = wspace.ene.epot_vdw;
	}

	void FF_SoftVDW::calcForces_NormalVDW()
	{
		if(!Active) return;

		// set up proxies to workspace to make code more readable.
		WorkSpace& wspace = getWSpace();
		size_t natom = wspace.atom.size();             // number of atoms in workspace
		const ParticleStore& atomparam = wspace.atom;  // atom parameter array
		SnapShotAtom *atom = wspace.cur.atom;						// atom coordinate array
		const NeighbourData *fnbor = wspace.nlist().getData(); // neighborlist

		int i, j, nj; // i&j are atom indices, nj count through neighbor list
		Tvector < double > fv; // force Tvector<double>
		double vdw_force = 0; // individual force magnitudes
		double vdw_potential = 0; // individual potential energy contributions

		double distij; // distance in Angstrom
		double invdistij, invd6;

		double sigma, epsilon, A, B;
		double qi;

		const double vdwinvSwidth = 1.0 / (VdwCutoff - VdwInnerCutoff);

		double vdwS;
		double vdwdSdd;

		double vdw14scaling;

		double atomi_radius;
		double atomi_epsilon;

		// reset the potential energy
		epot = 0;
		wspace.ene.epot_vdw = 0;
		wspace.ene.epot_vdw_att = 0;
		wspace.ene.epot_vdw_rep = 0;

		for(i = 0; i < natom; i++)
		{ 
			// loop over all particles
			atomi_radius = atomparam[i].radius;
			atomi_epsilon = atomparam[i].epsilon;
			qi = atomparam[i].charge;

			for(nj = 0; nj < fnbor[i].n; nj++) // and all their neighbours
			{
				j = NList32Bit_Index(fnbor[i].i[nj]);   // get atom index index
				if(j > i) break;                        // no shadow neighbors

				unsigned int bondOrder = NList32Bit_BondOrder(fnbor[i].i[nj]);
		
				if( bondOrder < BondOrder_1_4_Pair)
				{
					continue; // only get non-bonded and 1-4 neighbors 
				}
				else if( bondOrder == BondOrder_1_4_Pair )
				{
					vdw14scaling = Vdw14Scaling;
				}
				else
				{
					vdw14scaling = 1.0;
				}

				distij = wspace.atomxyz(i).dist( wspace.atomxyz(j) ); // in Angstrom, *NOT* SI units !!
				invdistij = 1.0 / distij;
				invd6 = invdistij * invdistij * invdistij * invdistij * invdistij * invdistij;

				// Van der Waals Force
				vdw_potential = 0.0;
				vdw_force = 0.0;

				if(distij < VdwCutoff)
				{
					// first calculate well depth (epsilon) and zero-force-separation (sigma)
					// Use Lorentz-Berthelot mixing rules
					// arthemic mean for sigma (=sum of radii) and geometric mean for epsilon
					sigma = atomi_radius + atomparam[j].radius;

					// two different mixing rules - at the moment only geometric mean is implemented
					epsilon = vdw14scaling * atomi_epsilon * atomparam[j].epsilon;

					B = sqr(sigma) * sqr(sigma) * sqr(sigma) * invd6;
					A = sqr(B); // within atom param !!!!

					vdw_potential = (2.0 * epsilon) * (0.5 * A - B);
					vdw_force = -(12.0 * epsilon * invdistij) * (A - B) / PhysicsConst::Angstrom;

					if(distij > VdwInnerCutoff)
					{
						vdwS = (1.0 - sqr(vdwinvSwidth * (distij - VdwInnerCutoff))); // Dimension less
						vdwdSdd = -4.0 * vdwinvSwidth * vdwinvSwidth * (distij - VdwInnerCutoff) * vdwS / PhysicsConst::Angstrom; // Units of per meter (not per angstrom!)
						vdwS *= vdwS;
						vdw_force = vdwS * vdw_force + vdw_potential * vdwdSdd;
						vdw_potential *= vdwS;
					}
				}

				// add up potentials -----------------------------

				wspace.ene.epot += vdw_potential;
				wspace.ene.epot_vdw += vdw_potential;
				if(vdw_potential > 0.0) wspace.ene.epot_vdw_rep += vdw_potential;
				else wspace.ene.epot_vdw_att += vdw_potential;

				// Now apply the force.
				fv.setTo(atom[i].p);
				fv.sub(atom[j].p); //fv is now the Tvector<double> along which all forces will act. fv points *from* atom i *to* atom nj
				fv.mul(invdistij * vdw_force);
				atom[i].f.sub(fv);
				atom[j].f.add(fv);
			}
		}

		epot = wspace.ene.epot_vdw;
	}
} // namespace Physics

