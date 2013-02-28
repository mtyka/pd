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
#include "protocolbase.h"
#include "maths/maths.h"
#include "workspace/workspace.h"
#include "forcefields/forcefield.h"
#include "workspace/neighbourlist.h"
#include "temperature.h"

namespace Protocol
{
	ProtocolBase::ProtocolBase(Physics::Forcefield & _ff )
		: WorkSpaceOperatorBase(_ff.getWSpace()),
		ff(&_ff),
		TargetTemp( new ConstantProfile( 300 ) ),
		OutputLevel(Verbosity::Normal) // All protocols talk normally by default
	{
		name = "Protocol";
		mon_counter = 0;
		nlist_counter = 0;
		UpdateNList = 1;
		settodefault();
	}

	void ProtocolBase::settodefault()
	{
		UpdateScr = 0;
		UpdateTra = 0;
		UpdateMon = 1;
		Steps = 0;
		mon_counter = 0;
	}

	int ProtocolBase::run()
	{
		if(Steps == 0)
		{
			throw(ArgumentException("You must set the number of steps before starting your simulation"));
		}
		if(ensureFFSetup()!=0) return -1;
		info();
		return runcore();
	}

	void ProtocolBase::listForcefields()
	{
		ff->listForcefields();
	}

	void ProtocolBase::printFinalStatistics(size_t steps)
	{
		printf("Timing: %ld sec (%d hrs %d mins %d secs)",
			(endtime - starttime),
			(endtime - starttime)/3600,
			((endtime - starttime)%3600)/60,
			((endtime - starttime)%3600)%60);
		if(steps>0)
		{
			printf("  %8.5f s/step", double(endtime - starttime)/double(steps));
		}
		printf("\n");
	}

	void ProtocolBase::addMonitor(Monitors::MonitorBase &monitor)
	{
		mon.add(monitor);
	}

	void ProtocolBase::clearMonitors()
	{
		mon.clear();
	}

	void ProtocolBase::popMonitor()
	{
		mon.pop_back();
	}

	void ProtocolBase::runmonitors()
	{
		if( UpdateMon == 0 ) return;
		if((mon_counter % UpdateMon)==0)
		{
			mon.measure();
		}
		mon_counter++;
	}

	void ProtocolBase::repeatmonitors()
	{
		if( UpdateMon == 0 ) return;
		if((mon_counter % UpdateMon)==0)
		{
			mon.repeat();
		}
		mon_counter++;
	}

	int ProtocolBase::ensureFFSetup()
	{
		ASSERT(ff!=NULL,CodeException,"ProtocolBase internal 'ff' pointer is invalid!");
		return ff->ensuresetup(getWSpace());
	}

    int ProtocolBase::ensureFFSetup(Physics::Forcefield& ff)
    {
        return ff.ensuresetup(getWSpace());
    }

	void ProtocolBase::assertStability()
	{
		if(!Maths::isNumber(getWSpace().ene.epot))
		{
			throw(ProcedureException("Simulation unstable or division by zero"));
		}
	}


	void ProtocolBase::refreshNeighborList()
	{
		if(UpdateNList<=0) return;
		if((nlist_counter%UpdateNList)==0)
		{
			getWSpace().cleanSpace();
			getWSpace().nlist().calcNewList();
		}
		nlist_counter++;
	}



	RangedProtocolBase::RangedProtocolBase( 
		Physics::Forcefield & _ff 
		): 
	ProtocolBase(_ff)
	{
		setFullRange();
	}

	RangedProtocolBase::RangedProtocolBase( 
		Physics::Forcefield & _ff, 
		const PickAtomRange& _newRange 
		):
	ProtocolBase(_ff)
	{
		setRange(_newRange);
	}

	void RangedProtocolBase::setFullRange() 
	{ 
		m_Range = PickAtomRange( getWSpace() ); 
	}


	RangesProtocolBase::RangesProtocolBase( 
		Physics::Forcefield & _ff
		)
		: ProtocolBase(_ff),
		m_PickChangeSerial(0)
	{
		setFullRange();
	}

	RangesProtocolBase::RangesProtocolBase( 
		Physics::Forcefield & _ff, 
		const PickAtomRange& _newRange
		)
		: ProtocolBase(_ff),
		m_PickChangeSerial(0)
	{
		setRange(_newRange);
	}

	RangesProtocolBase::RangesProtocolBase( 
		Physics::Forcefield& _ff, 
		const PickAtomRanges& _newRange
		)
		: ProtocolBase(_ff),
		m_PickChangeSerial(0)
	{
		setRange(_newRange);
	}

	void RangesProtocolBase::setRange(const PickAtomRanges& _newRange)
	{	
		m_PickChangeSerial++;
		m_Range.clear();
		m_Range = _newRange;
	}

	void RangesProtocolBase::setFullRange() 
	{ 
		m_PickChangeSerial++;
		m_Range.clear();
		m_Range.addRange(PickAtomRange( getWSpace() ) ); 
	}


	PickedProtocolBase::PickedProtocolBase( 
		Physics::Forcefield& _ff
		) 
		: ProtocolBase(_ff),
		m_PickerSerial(0)
	{
		setPickEverything();
	}

	PickedProtocolBase::PickedProtocolBase( 
		Physics::Forcefield& _ff, 
		const PickBase& _Picker
		) 
		: ProtocolBase(_ff),
		m_PickerSerial(0)
	{
		setPicking(_Picker);
	}

	void PickedProtocolBase::setPicking(const PickBase& _Picker)
	{
		m_PickerSerial++;
		m_Atoms.setPicking( getWSpace(),_Picker);
		m_Picker = _Picker;
	}

	void PickedProtocolBase::setPickEverything()
	{
		m_PickerSerial++;
		PickAllParticles picker;
		m_Atoms.setPicking( getWSpace(),picker);
		m_Picker = picker;
	}



	bool every(int _Step, int _Interval)
	{
		if(_Interval <= 0) return false;
		if(_Step <= 0)     return false;
		if(((_Step) % _Interval) == 0) return true;
		return false;
	}




	ProtocolSet::ProtocolSet( Physics::Forcefield & _ff ):
		ProtocolBase( _ff )
	{
		name = "ProtocolSet[empty]";	
	}

	int ProtocolSet::run()
	{
		int result = 0;
		for(size_t i = 0; i < size(); i++) result += element(i).run();
		return result;
	}


	int ProtocolSet::runcore()
	{
		int result = 0;
		for(size_t i = 0; i < size(); i++) result += element(i).runcore();
		return result;
	}

	/// prints a little block of parameter information	
	void ProtocolSet::info() const
	{
		for(size_t i = 0; i < size(); i++) element(i).info();
	}

	void ProtocolSet::printFinalStatistics(size_t steps)
	{
		for(size_t i = 0; i < size(); i++) element(i).printFinalStatistics( steps );
	}

	/// Add monitors to the simulation, just pass a reference to any monitor
	void ProtocolSet::addMonitor(Monitors::MonitorBase &monitor)
	{
		for(size_t i = 0; i < size(); i++) element(i).addMonitor( monitor );
	}

	/// Remove all monitors from this simulation
	void ProtocolSet::clearMonitors()
	{
		for(size_t i = 0; i < size(); i++) element(i).clearMonitors(  );
	}

	/// remove the last monitor
	void ProtocolSet::popMonitor()
	{
		for(size_t i = 0; i < size(); i++) element(i).popMonitor(  );
	}

	/// Temperature (only if Thermostats are active) - note this does not take a number -
	/// it takes a Temperature type, i.e. Temp, TempLinear or TempExp. These allow the
	/// user to define different temperature protocols.

	void ProtocolSet::setTargetTemp( const ProfileBase & _TargetTemp ) 
	{
		for(size_t i = 0; i < size(); i++) element(i).setTargetTemp( _TargetTemp );
	}

	void ProtocolSet::setTargetTemp( double _Temperature )
	{
		for(size_t i = 0; i < size(); i++) element(i).setTargetTemp( _Temperature );
	}

	void ProtocolSet::generateName(){
		
		name = "ProtocolSet[";
		if( size() == 0 ){ 
			name += "emtpy";
		}else{
			for(size_t i = 0; i < size(); i++)
			{
				name += element(i).name;
				if( i < (size()-1)) name += ",";
			}
		}
		name += "]";
	}


}














