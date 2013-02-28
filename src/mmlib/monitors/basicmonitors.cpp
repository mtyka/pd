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

#include "maths/maths_vector.h"
#include "workspace/workspace.h"
#include "workspace/space.h"
#include "workspace/neighbourlist.h"
#include "forcefields/forcefield.h"
#include "forcefields/ffbonded.h"   // Provides class member 'dehedral'
#include "forcefields/restraintbase.h"  
#include "filters/filterbase.h" 

#include "monitors/basicmonitors.h"

namespace Monitors
{
	PotentialEnergyMonitor::PotentialEnergyMonitor(const WorkSpace& _wspace )
	{
		name = "monepot";
		wspace = &_wspace;
	}

	void PotentialEnergyMonitor::setcurdata()
	{
		addData(  wspace->ene.epot * Physics::PhysicsConst::J2kcal * Physics::PhysicsConst::Na );
	}


	VolumeMonitor::VolumeMonitor(const ClosedSpace& _space )
	{
		name = "monvolume";
		space = &_space;
	}

	void VolumeMonitor::setcurdata()
	{
		addData(  space->volume() ); 
	}


	CRMSMonitor::CRMSMonitor(const WorkSpace& _wspace)
	{
		name = "crms";
		wspace = &_wspace;
	}

	void CRMSMonitor::setcurdata()
	{
		addData(  wspace->calcCRMS_CA() );
	}


	DRMSMonitor::DRMSMonitor(const WorkSpace& _wspace)
	{
		name = "drms";
		wspace = &_wspace;
	}

	void DRMSMonitor::setcurdata()
	{
		THROW(CodeException,"Temporarily broken code\n");
		//addData(  wspace->calcdRMS() );
	}


	DihedralMonitor::DihedralMonitor(const WorkSpace& _wspace, const IndexQuartet &_dihedral)
	{
		name = "Dihed Mon";
		wspace = &_wspace;
		dihedral = _dihedral;

		if((dihedral.i<0)||dihedral.i>=wspace->atom.size()){
			printf("dihedral.i = %d (wspace->atom.size() = %d ) \n",dihedral.i,wspace->atom.size());
			THROW(ArgumentException,"Dihedral atom i is out of range for particle system");
		}
		if((dihedral.a<0)||dihedral.a>=wspace->atom.size()){
			printf("dihedral.a = %d (wspace->atom.size() = %d ) \n",dihedral.a,wspace->atom.size());
			THROW(ArgumentException,"Dihedral atom a is out of range for particle system");
		}
		if((dihedral.b<0)||dihedral.b>=wspace->atom.size()){
			printf("dihedral.b = %d (wspace->atom.size() = %d ) \n",dihedral.b,wspace->atom.size());
			THROW(ArgumentException,"Dihedral atom b is out of range for particle system");
		}
		if((dihedral.j<0)||dihedral.j>=wspace->atom.size()){
			printf("dihedral.j = %d (wspace->atom.size() = %d ) \n",dihedral.j,wspace->atom.size());
			THROW(ArgumentException,"Dihedral atom j is out of range for particle system");
		}
	}

	void DihedralMonitor::setcurdata(){
		curdata =
			Maths::calcTorsionAngle(
				wspace->cur.atom[dihedral.i].p,
				wspace->cur.atom[dihedral.a].p,
				wspace->cur.atom[dihedral.b].p,
				wspace->cur.atom[dihedral.j].p );
	}


	ForcefieldEnergyMonitor::ForcefieldEnergyMonitor(const Physics::ForcefieldBase& newff )
	{
		name = "monff";
		ff = &newff;
	}

	void ForcefieldEnergyMonitor::setcurdata()
	{
		addData(  ff->epot * Physics::PhysicsConst::J2kcal * Physics::PhysicsConst::Na );
	}


	RestraintMonitor::RestraintMonitor(const Physics::RestraintForcefieldBase& newff)
	{
		name = "monrest";
		ff = &newff;
	}

	void RestraintMonitor::setcurdata()
	{
		addData(  ff->deviat );
	}



	AtomDistanceMonitor::AtomDistanceMonitor(const MoleculeBase& _wspace, unsigned _atomi, unsigned _atomj ){
		name = "mondist";
		wspace = &_wspace;
		Atom_i = _atomi;
		Atom_j = _atomj;
		if(Atom_i>=wspace->nAtoms()) throw(ArgumentException("Atom Index i out of range"));
		if(Atom_j>=wspace->nAtoms()) throw(ArgumentException("Atom Index j out of range"));
	}

	void AtomDistanceMonitor::setcurdata()
	{
		addData(  Maths::dist(wspace->atomxyz(Atom_i),wspace->atomxyz(Atom_j)) );
	}



	DistributionFunctionMonitor::DistributionFunctionMonitor(
		const WorkSpace& _wspace, 
		const PickBase &picker1, 
		const PickBase &picker2
	){
		name = "mondistfunc";
		wspace = &_wspace;
		m_Cutoff = 10;
		//m_Update = 10;
		setSet1(picker1);
		setSet2(picker2);

	}

	void DistributionFunctionMonitor::setSet1( const PickBase &picker){
		// set the atom selection
		m_Set1 = SelectedIndices( *wspace, picker);
	}


	void DistributionFunctionMonitor::setSet2( const PickBase &picker){
		// set the atom selection
		m_Set2 = SelectedIndices( *wspace, picker);
	}


	void DistributionFunctionMonitor::setcurdata()
	{

		
		const NeighbourListBase &nlist = wspace->nlist();

		SelectedIndices *seta = &m_Set2;
		SelectedIndices *setb = &m_Set1;

		for(size_t loop=0; loop < 2; loop++ )
		{
			if( loop == 1){
				seta = &m_Set1;
        setb = &m_Set2;
			}
			for(size_t si=0; si < seta->size(); si++ )
			{
				size_t i = (*seta)[si];
				// use neighbourlist to loop over neighbours
				for(int nj = 0; nj < nlist.nNeighbours(i); nj++)
				{	
			
					int j = nlist.getNeighbourIndex(i,nj);

					if(j>=i) continue;

					// this could be a better search provided by Selected Indices
					size_t tj;
					for(tj=0; tj < setb->size(); tj++ ){ 
						if(j== (*setb)[tj]) break;				
					}
					if( tj>=setb->size() ) continue; // does not match m_Set2;

					Maths::dvector dc;
					dc.diff( wspace->cur.atom[i].p, wspace->cur.atom[j].p );

					wspace->boundary().getClosestImage(dc);

					double sqrdist = dc.innerdot();
					if(sqrdist < Maths::sqr(m_Cutoff)){
						double distij = sqrt(sqrdist);
						addData(distij);
					}
				}

			}

		}
		
	}






	TIMonitor::TIMonitor(const Physics::FF_Extension_TI& newti )
	{
		name = "monti";
		ti = &newti;
	}

	void TIMonitor::setcurdata()
	{
		addData(  ti->get_dEdlambda() );
	}


	FilterMonitor::FilterMonitor()
	{
		m_Filter = NULL;
	}

	FilterMonitor::FilterMonitor(FilterBase &_Filter)
	{
		setFilter(_Filter);
	}

	void FilterMonitor::setFilter(FilterBase &_Filter)
	{
		m_Filter = &_Filter;
	}
	
	void  FilterMonitor::setcurdata()
	{
		addData(  m_Filter != NULL && m_Filter->passes() ? 1.0 : 0.0 );
	}

} // Namespace Monitors

