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

#ifndef __MONITORS_H
#define __MONITORS_H

// Essential Headers
#include "monitors/monitorbase.h" // Provides a base class
#include "workspace/workspace.fwd.h"

// Forward declarations
namespace Physics
{
	class PD_API ForcefieldBase;	
	class PD_API RestraintForcefieldBase;
	class PD_API FF_Extension_TI;
}

class PD_API FilterBase;
class PD_API MoleculeBase;
class PD_API ClosedSpace;
class PD_API PickBase;

namespace Monitors
{
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
	class PD_API PotentialEnergyMonitor : public MonitorBase
	{
	public:
		PotentialEnergyMonitor(const WorkSpace& _wspace );
		virtual PotentialEnergyMonitor* clone() const { return new PotentialEnergyMonitor(*this); }

	protected:
		void setcurdata();

	private:
		const WorkSpace *wspace;
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
	class PD_API VolumeMonitor : public MonitorBase
	{
	public:
		VolumeMonitor(const ClosedSpace& _space);
		virtual VolumeMonitor* clone() const { return new VolumeMonitor(*this); }

	protected:
		void setcurdata();

	private:
		const ClosedSpace *space;
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
	class PD_API ForcefieldEnergyMonitor: public MonitorBase
	{
	public:
		ForcefieldEnergyMonitor(const Physics::ForcefieldBase& newff);
		virtual ForcefieldEnergyMonitor* clone() const { return new ForcefieldEnergyMonitor(*this); }

	protected:
		void setcurdata();

	private:
		const Physics::ForcefieldBase *ff;
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
	class PD_API RestraintMonitor: public MonitorBase
	{
	public:
		RestraintMonitor(const Physics::RestraintForcefieldBase& newff );
		virtual RestraintMonitor* clone() const { return new RestraintMonitor(*this); }
		void setcurdata();

	private:
		const Physics::RestraintForcefieldBase *ff;
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
	class PD_API AtomDistanceMonitor: public MonitorBase
	{
	public:
		AtomDistanceMonitor(const MoleculeBase& _wspace, unsigned _atomi, unsigned _atomj );
		virtual AtomDistanceMonitor* clone() const { return new AtomDistanceMonitor(*this); }
		void setcurdata();

	private:
		unsigned Atom_i;
		unsigned Atom_j;
		const MoleculeBase *wspace;
	};


	//-------------------------------------------------
	//
	/// \brief  Class to measure distances between sets of atoms. This is useful to measure
	///         distribution functions for example.
	///
	/// \details Will measure distances between all the atoms defined in set 1 vs all the atoms
	///         defined in set two, provided atom index i<j. (this avoids double counting and
	///         delf-distances. The two sets are not mutually exclusive which makes this Monitor
	///         very powerful when combined with interesting PickBases. Note that for every distance
	///         a Sqrt must be evaluated which is slow! THus this monitor should be used with care
	///         or post-simulation to avoid slowing the simulation excessively.
	///         A distance limit can be set to include only close interactions ! (Cutoff)
	///         Cutoff is set to 6A by default.
	///
	/// \author Mike Tyka
	///
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class PD_API DistributionFunctionMonitor: public MonitorBase
	{
	public:
		DistributionFunctionMonitor(const WorkSpace& _wspace, const PickBase &picker1, const PickBase &picker2);
		virtual DistributionFunctionMonitor* clone() const { return new DistributionFunctionMonitor(*this); }
		void setcurdata();

		void setSet1( const PickBase &picker);	
		void setSet2( const PickBase &picker);	

		void   setCutoff(double _cutoff){ m_Cutoff = _cutoff; }
		double getCutoff() const { return m_Cutoff; }

	private:	
		SelectedIndices m_Set1;
		SelectedIndices m_Set2;

		double m_Cutoff;
		const WorkSpace *wspace;
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
	class PD_API TIMonitor: public MonitorBase
	{
	public:
		TIMonitor(const Physics::FF_Extension_TI& newti );
		virtual TIMonitor* clone() const { return new TIMonitor(*this); }
		void setcurdata();

	private:
		const Physics::FF_Extension_TI *ti;
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
	class PD_API CRMSMonitor: public MonitorBase
	{
	public:
		CRMSMonitor(const WorkSpace& _wspace);
		virtual CRMSMonitor* clone() const { return new CRMSMonitor(*this); }
		void setcurdata();

	private:
		const WorkSpace* wspace;
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
	class PD_API DRMSMonitor: public MonitorBase{
	public:

		DRMSMonitor(const WorkSpace& _wspace);
		virtual DRMSMonitor* clone() const { return new DRMSMonitor(*this); }

		void setcurdata();
	private:
		const WorkSpace* wspace;
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
	class PD_API DihedralMonitor: public MonitorBase
	{
	public:
		DihedralMonitor(const WorkSpace& _wspace, const IndexQuartet &_dihedral);
		virtual DihedralMonitor* clone() const { return new DihedralMonitor(*this); }
		void setcurdata();

	private:
		const WorkSpace *wspace;
		IndexQuartet dihedral;
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
	class PD_API FilterMonitor : public MonitorBase
	{
	public:
		FilterMonitor();
		FilterMonitor(FilterBase &_Filter);
		virtual FilterMonitor* clone() const { return new FilterMonitor(*this); }
		void setFilter(FilterBase &_Filter);

	protected:
		virtual void  setcurdata();
		FilterBase *m_Filter;
	};
}

#endif
