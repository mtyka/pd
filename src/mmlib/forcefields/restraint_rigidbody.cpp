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

#include "forcefields/restraint_rigidbody.h"
#include "workspace/neighbourlist.h"
#include "workspace/space.h"

using namespace Maths;

namespace Physics
{



// -------------------------------------------------------------------------
//
//    Translational restraint
//

double FF_BodyRestraint_Position::calcEnergyAtQ(double newQ)
{ 
	return 0.0; 
}

void FF_BodyRestraint_Position::info() const
{
}

void FF_BodyRestraint_Position::setup()
{
	FF_BodyRestraint::setup();
}


void  FF_BodyRestraint_Position::calcForces()
{
	// Proxies
	WorkSpace& wspace = getWSpace();

	size_t i;
	double ene;
	dvector force;

	double k_SI = k / (PhysicsConst::J2kcal * PhysicsConst::Na); // convert to SI units for calculation
	double forcemag;
	double dist;

	epot = 0;

	// first calculate the current centre of mass/gravity (depending on weight setting)
	// and the restraint position one
	dvector c_centre(0,0,0);
	dvector r_centre(0,0,0);
	double  totalweight = 0;
	dvector temp;	
	for(i = 0; i < m_Selection.size(); i++) 
	{
		// get weight
		double imass = 1.0;
		if(MassWeight) imass = wspace.atom[m_Selection.iat(i)].mass;

		//printf("%e \n",imass);
		// save positions
		temp.setTo(m_Selection.p(i));
		temp.mul(imass);
		c_centre.add(temp);

		temp.setTo(m_Selection.sp(i));
		temp.mul(imass);
		r_centre.add(temp);

		totalweight += imass;
	}

	c_centre.div(totalweight);
	r_centre.div(totalweight);

	// "force" on centre of mass is calculated 
	force.diff(c_centre, r_centre);
	wspace.boundary().getClosestImage(force);
	dist = force.mag();

	// multiply by distance appropritae number of times for the power of 
	// the restraint (default=2)
	forcemag = 0.5 * k_SI;
	for(int p = 0; p < (Power - 2); p++) {
		forcemag *= dist;
	}

	ene = forcemag;
	ene *= dist * dist;
	forcemag *= double (Power);
	force.mul(-forcemag / PhysicsConst::Angstrom);
	epot += ene;
	
	dvector iforce;
	double invtotalmass = 1.0/totalweight;

	if(!Passive) {
		wspace.ene.epot += ene;

		for(i = 0; i < m_Selection.size(); i++) {
			// get weight
			double imass = 1.0;
			if(MassWeight) imass = wspace.atom[m_Selection.iat(i)].mass;
	
			m_Selection.f(i).x += force.x * invtotalmass * imass;
			m_Selection.f(i).y += force.y * invtotalmass * imass;
			m_Selection.f(i).z += force.z * invtotalmass * imass;
		}

	}
	deviat = 2.0*epot/k_SI;
}



// -------------------------------------------------------------------------
//
//    Body Angle Restraint
//

double FF_BodyRestraint_Angle::calcEnergyAtQ(double newQ)
{ 
	return 0.0; 
}

void FF_BodyRestraint_Angle::info() const
{
	printf("Force Constant:     %8.4f\n",k);
	printf("Form:      ");
	if(Form == Harmonic) printf("Harmonic\n");
	if(Form == Cosine) printf("Cosine\n");
}

void FF_BodyRestraint_Angle::setup()
{
	FF_BodyRestraint::setup();
	WorkSpace& wspace = getWSpace();

	// calculate three centre of masses of the target structure

	dvector centre[3] = { dvector(0,0,0), dvector(0,0,0), dvector(0,0,0) };
	dvector fcentre[3];
	double  totalweight[3] = {0,0,0};
	size_t  firstthird = m_Selection.size()/3;
	size_t  secondthird = firstthird*2;
	size_t  istart[3] = {0,firstthird,secondthird};
	size_t  iend[3]   = {firstthird+1,secondthird+1,m_Selection.size()};
	size_t  ithird;
	printf("Setting Restraint Pos\n");
	dvector temp;	
	int i;
	for(ithird=0;ithird<3;ithird++)
	{
		for(i = istart[ithird]; i < iend[ithird]; i++) 
		{
			double imass = 1.0;
			if(MassWeight) imass = wspace.atom[m_Selection.iat(i)].mass;
			temp.setTo(m_Selection.sp(i));
			temp.mul(imass);
			centre[ithird].add(temp);
			totalweight[ithird] += imass;
		}
		centre[ithird].div(totalweight[ithird]);
	}

	m_AngleTerm_Vector.setTo(centre[0]);
	m_AngleTerm_Vector.sub(centre[1]);
	m_AngleTerm_Vector.unify();

}

void  FF_BodyRestraint_Angle::calcForces()
{
	// Proxies
	WorkSpace& wspace = getWSpace();

	size_t i;
	double ene;
	dvector force;

	double k_SI = k / (PhysicsConst::J2kcal * PhysicsConst::Na); // convert to SI units for calculation

	epot = 0;

	// first calculate three centre of masses

	dvector centre[3] = { dvector(0,0,0), dvector(0,0,0), dvector(0,0,0) };
	dvector fcentre[3];
	double  totalweight[3] = {0,0,0};
	size_t  firstthird = m_Selection.size()/3;
	size_t  secondthird = firstthird*2;
	size_t  istart[3] = {0,firstthird,secondthird};
	size_t  iend[3]   = {firstthird+1,secondthird+1,m_Selection.size()};
	size_t  ithird;

	dvector temp;	
	for(ithird=0;ithird<3;ithird++){
		for(i = istart[ithird]; i < iend[ithird]; i++) {
			double imass = 1.0;
			if(MassWeight) imass = wspace.atom[m_Selection.iat(i)].mass;
			temp.setTo(m_Selection.p(i));
			temp.mul(imass);
			centre[ithird].add(temp);
			totalweight[ithird] += imass;
		}
		centre[ithird].div(totalweight[ithird]);
	}


	// Angle term
	// Restraint on angle between the vector between the first and second centre of mass
	// and an arbitrary axis
	dvector iv, jv, uiv, ujv;   // bond vectors

	double force_magnitude;
	dvector forcev;             // force dvector 
	dvector bondnormal;         // normal dvector  (plane of bond)

	iv.setTo(centre[0]);
	iv.sub(centre[1]);	
	jv.setTo(m_AngleTerm_Vector); // unit vector
	
	double a = 1.0 / (iv.mag()*jv.mag());
  double t = iv.x*jv.x + iv.y*jv.y + iv.z*jv.z;
	double P = a*t;
	if(P>1.00) P=1.0;
	if(P<-1.00) P=-1.0;
	double phi = acos(P);
	ene = 0;
	force_magnitude = 0;
	double E=0;
	double dEdP=0;
	if(Form == Harmonic){
		E = 0.5*k_SI*sqr(phi);
		if(P==1.0) dEdP = 0;
		else       dEdP = k_SI*phi/sqrt(1.0 - sqr(P))/ (PhysicsConst::Angstrom);
	} else
	if(Form == Cosine){
		E = 0.5 * k_SI * (1.0-P);
	  dEdP = 0.5 * k_SI / (PhysicsConst::Angstrom);
	}
	
	dvector dPdi(iv);
	dPdi.mul(-t/iv.innerdot());
	dPdi.add(jv);
	dPdi.mul(a);

	fcentre[0].setTo(dPdi);
	fcentre[0].mul(dEdP);
	fcentre[1].setTo(fcentre[0]);
	fcentre[1].mul(-1);

	epot += E;
	if(!Passive) {
		wspace.ene.epot += E;
	
		for(ithird=0;ithird<2;ithird++){
			for(i = istart[ithird]; i < iend[ithird]; i++) {
				// get atom weight
				double imass = 1.0;
				if(MassWeight) imass = wspace.atom[m_Selection.iat(i)].mass;
	
				m_Selection.f(i).x += fcentre[ithird].x * imass / totalweight[ithird];
				m_Selection.f(i).y += fcentre[ithird].y * imass / totalweight[ithird];
				m_Selection.f(i).z += fcentre[ithird].z * imass / totalweight[ithird];
			}
		}

	}
	deviat = 2.0*epot/k_SI;
	Q = phi;
}












double FF_BodyRestraint_Dihedral::calcEnergyAtQ(double newQ)
{ 
	return 0.0; 
}

void FF_BodyRestraint_Dihedral::info() const
{
	printf("Force Constant:     %8.4f\n",k);
	printf("Form:      ");
	if(Form == Harmonic) printf("Harmonic\n");
	if(Form == Cosine) printf("Cosine\n");
}

void FF_BodyRestraint_Dihedral::setup()
{
	FF_BodyRestraint::setup();
	WorkSpace& wspace = getWSpace();

// calculate three centre of masses of the target structure

	dvector centre[3] = { dvector(0,0,0), dvector(0,0,0), dvector(0,0,0) };
	dvector fcentre[3];
	double  totalweight[3] = {0,0,0};
	size_t  firstthird = m_Selection.size()/3;
	size_t  secondthird = firstthird*2;
	size_t  istart[3] = {0,firstthird,secondthird};
	size_t  iend[3]   = {firstthird+1,secondthird+1,m_Selection.size()};
	size_t  ithird;

	dvector temp;	
	int i;
	for(ithird=0;ithird<3;ithird++){
		for(i = istart[ithird]; i < iend[ithird]; i++) {
			double imass = 1.0;
			if(MassWeight) imass = wspace.atom[m_Selection.iat(i)].mass;
			temp.setTo(m_Selection.sp(i));
			temp.mul(imass);
			centre[ithird].add(temp);
			totalweight[ithird] += imass;
		}
		centre[ithird].div(totalweight[ithird]);
	}

	m_DihedralTerm_Vector.diff(centre[2],centre[1]);
	m_DihedralTerm_Vector.unify();

}

void  FF_BodyRestraint_Dihedral::calcForces()
{
	// Proxies
	WorkSpace& wspace = getWSpace();

	size_t i;
	dvector force;

	double k_SI = k / (PhysicsConst::J2kcal * PhysicsConst::Na); // convert to SI units for calculation

	epot = 0;

	// first calculate three centre of masses
	dvector centre[3] = { dvector(0,0,0), dvector(0,0,0), dvector(0,0,0) };
	dvector fcentre[3];
	double  totalweight[3] = {0,0,0};
	size_t  firstthird = m_Selection.size()/3;
	size_t  secondthird = firstthird*2;
	size_t  istart[3] = {0,firstthird,secondthird};
	size_t  iend[3]   = {firstthird+1,secondthird+1,m_Selection.size()};
	size_t  ithird;

	dvector temp;	
	for(ithird=0;ithird<3;ithird++){
		for(i = istart[ithird]; i < iend[ithird]; i++) {
			double imass = 1.0;
			if(MassWeight) imass = wspace.atom[m_Selection.iat(i)].mass;
			temp.setTo(m_Selection.p(i));
			temp.mul(imass);
			centre[ithird].add(temp);
			totalweight[ithird] += imass;
		}
		centre[ithird].div(totalweight[ithird]);
	}

	// Dihedral Term

	dvector e,k,g,h,exk,gxk;
	e.setTo(m_DihedralTerm_Vector);
	k.diff(centre[0],centre[1]);
	g.diff(centre[2],centre[1]);
	h.diff(centre[0],centre[2]);
	exk.crossProduct(e,k);
	gxk.crossProduct(g,k);

	double a = 1.0 / (exk.mag()*gxk.mag());
  double t = exk.x*gxk.x + exk.y*gxk.y + exk.z*gxk.z;
	double P = a*t;
	if(P>1.00) P=1.0;
	if(P<-1.00) P=-1.0;
	double phi = acos(P);
	//printf("Qd: %f \n",phi);
	double E = 0.5*k_SI*sqr(phi);
	double dEdP;
	if(P==1.0) dEdP = 0;
	else       dEdP = k_SI*phi/sqrt(1.0 - sqr(P));

	//printf("%e %e %e %e %e %e\n",a,t,P,phi,E,dEdP);
	dvector dPdexk(exk);
	dPdexk.mul(-t/exk.innerdot());
	dPdexk.add(gxk);
	dPdexk.mul(a);

	dvector dPdgxk(gxk);
	dPdgxk.mul(-t/gxk.innerdot());
	dPdgxk.add(exk);
	dPdgxk.mul(a);

	dvector part2;
	dvector dPda;
	dPda.crossProduct(dPdexk,e);
	part2.crossProduct(dPdgxk,g);
	dPda.add(part2);

	dvector dPdb;
	dPdb.crossProduct(e,dPdexk);
	part2.crossProduct(dPdgxk,h);
	dPdb.add(part2);
	
	dvector dPdc;
	dPdc.crossProduct(dPdgxk,k);
	dPdc.mul(-1);

	dvector dEdc[3];
	dEdc[0].setTo(dPda); dEdc[0].mul(dEdP);
	dEdc[1].setTo(dPdb); dEdc[1].mul(dEdP);
	dEdc[2].setTo(dPdc); dEdc[2].mul(dEdP);

	
	epot += E;
	if(!Passive) {
		wspace.ene.epot += E;
	
		for(ithird=0;ithird<3;ithird++){
			for(i = istart[ithird]; i < iend[ithird]; i++) {
				// get atom weight
				double imass = 1.0;
				if(MassWeight) imass = wspace.atom[m_Selection.iat(i)].mass;
	
				m_Selection.f(i).x += dEdc[ithird].x * imass / (totalweight[ithird]*PhysicsConst::Angstrom);
				m_Selection.f(i).y += dEdc[ithird].y * imass / (totalweight[ithird]*PhysicsConst::Angstrom);
				m_Selection.f(i).z += dEdc[ithird].z * imass / (totalweight[ithird]*PhysicsConst::Angstrom);
			}
		}
	}

	deviat = 2.0*epot/k_SI;
	Q = phi;
}


}

