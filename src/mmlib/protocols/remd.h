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

#ifndef __REMD_H
#define __REMD_H

// Essential Headers
#include <vector>
#include "workspace/snapshot.h" // Provides a class member
#include "workspace/workspace.fwd.h"


namespace Protocol{
    class PD_API ProtocolBase;
    class PD_API REX_Local;






    //-------------------------------------------------
    //
    /// \brief  Holds a REX_Replica for REX_ReplicaExchange Simulations 
    ///
    /// \details 
    ///
    /// \author Mike Tyka 
    ///
    ///
    class PD_API REX_Replica{
        friend class PD_API REX_Local;

    public:	
        enum ExchangeModeType { Swap, Downward, Upward };

        REX_Replica(const Protocol::ProtocolBase & _protocol);

        /// Sets the temperature at which this replica should run. 
        void setTargetTemp( float _Temperature ){
            m_Temperature = _Temperature;
            protocol().setTargetTemp( m_Temperature );
        }

        /// Non-Constant Accessor for the protocol handling this replica
        Protocol::ProtocolBase & protocol()       { return m_Protocol.data(); };

        /// Constant Accessor for the protocol handling this replica
        const Protocol::ProtocolBase & protocol() const { return m_Protocol.data(); };

        /// Returns the target temperature of this REX_Replica.
        double getTargetTemp() const { return m_Temperature; };

        /// Sets the velocities of the internal SnapShot to be a Maxwell distribution at
        /// a temperature of m_Temperature. 
        void setInitialSpeeds();

        /// Print a brief description of this object
        void info() const {
            printf("Rep   Sim=%s  T= %5.1f  ExchangeMode=", protocol().name.c_str(), getTargetTemp() );
            switch ( ExchangeMode ){
                case REX_Replica::Swap:
                    printf("Swap");
                    break;
                case REX_Replica::Downward:
                    printf("Downward");
                    break;
                case REX_Replica::Upward:
                    printf("Upward");
                    break;
            }
            printf("\n");
        }

        /// This controls how this replica exchanges with the replica <i> above </i>. Possible values are Swap (Swaps the structures),
        /// Downward (This replica is set to a copy of the one above, the upper one remains unchanged) and Upward (The upper replica
        /// is set to this one, while thisone remains unchanged).
        ExchangeModeType ExchangeMode;

    private:

        /// The target temperature of this replica
        double m_Temperature;

        /// The simulation/algorithm currently responsible for simulating this replica.

        /// This holds the simulation protocol which simulates replicas at this temperature.
        CloneHolder< Protocol::ProtocolBase > m_Protocol;

        /// The current state of the system
        SnapShot m_Structure;
    };








    //-------------------------------------------------
    //
    /// \brief Local REX_Replica Exchange Dynamics
    ///
    /// \details 
    ///    
    /// REX_ReplicaExchange Molecular Dynamics:
    /// Yuji Sugita, Yuko Okamoto, REX_Replica-exchange molecular dynamics method
    /// for protein folding. Chemical Physics Letters 314 141-151 (1999)
    ///
    /// \author Mike Tyka 
    ///
    class PD_API REX_Local:public ProtocolBase
    {
    public:
        REX_Local(Physics::Forcefield & _ff):
          ProtocolBase(_ff)
          {
              FocusRep = 0;
              name = "Local Replica Exchange Manager";
          }

          virtual ~REX_Local(){
          }

          virtual REX_Local* clone() const 
          { 
              return new REX_Local(*this); 
          }

          virtual int runcore();

          /// prints a little block of parameter information
          virtual void info() const{ 
              printf("---- %s ------------------\n", name.c_str());
              printf("Steps(Rounds)   %d \n", Steps );
              printf("UpdateScr       %d \n", UpdateScr );
              printf("UpdateTra       %d \n", UpdateTra );
              printf("UpdateMon       %d \n", UpdateMon );
              printf("FocusRep        %d \n", FocusRep );
              printf("Replicas:       %d \n", rep.size() );
              for(size_t r=0;r<rep.size();r++){
                  printf("Rep  %d : ", r );
                  rep[r].info();
              }
          }

          /// prints a little block of parameter information
          virtual void detail() const{ 
              info();
              for(size_t r=0;r<rep.size();r++){
                  rep[r].protocol().info();	
              }
          }


          /// Grabs a copy of all the structures 
          void getReplicas(SnapShot * m_Structure_ext);

          /// Sets the structure of every replica to the snapshot provided
          void setAllReplicasTo(SnapShot &m_Structure_ext);

          /// Sets the structures of all the replicas to those in the array of snapshots provided
          void setReplicas(SnapShot * m_Structure_ext);

          /// \brief Adds a replica by using the passed simulation as a template
          /// Note that an internal copy is made, i.e. no references to the template
          /// are retained.
          void addReplica(
              const Protocol::ProtocolBase& simTemplate, 
              float Temperature,
              REX_Replica::ExchangeModeType exmode = REX_Replica::Swap
              );

          /// \brief Adds replicas by using the passed MD simulation as a template
          /// Note that an internal copy is made, i.e. no references to the template
          /// are retained.
          void addReplicas(
              const Protocol::ProtocolBase &simTemplate, 
              std::vector<float> newtemps,
              REX_Replica::ExchangeModeType exmode = REX_Replica::Swap
              );

          /// \brief Adds replicas by using the passed MD simulation as a template
          /// Note that an internal copy is made, i.e. no references to the template
          /// are retained. This version adds N replicas at exponentially increasing temp (each time multiplying the
          /// previous temperature by 'factor', starting with firstTmp
          void addReplicas(
              const Protocol::ProtocolBase &simTemplate, 
              float firstTemp,
              float factor,
              size_t N,
              REX_Replica::ExchangeModeType exmode = REX_Replica::Swap
              );

          /// Removes all the added replicas
          void clearReplicas()
          {
              rep.clear();
          };

          /// prints a MIME encoded ascii print of all the replcas in order 
          void printCheckPointMIME();

          /// this function can be used to restart a simulation given a file that
          /// containes all the replica snapshots.
          void readCheckPointMIME(const std::string &filename);

    public: 

        /// This variable determines which replica the Monitors added to this protocol
        /// are watching. It is simply the number of the replica, starting with 0.
        int FocusRep;

    protected:

        /// This function is called to decde if two replias should be exchanged in some way
        /// This function can be virtually overloaded to create algorithms with different criteria
        /// but by default this uses the metropolis criterion from the REMD algorithm introduced by Sugita et al.
        virtual bool exchangeCriterion( const REX_Replica &lowerrep,
            const REX_Replica &upperrep ) const;

        /// CheckInternal state - 
        /// Check that all replicas have identical numbers of atoms and work with the same workspace
        void doInternalSafetyCheck();

        /// prints a line of current energies
        virtual void infoLine() const {}; 

        /// prints the headers for the above function
        virtual void infoLineHeader() const{}; 

        std::vector <REX_Replica> rep;
    };
}
#endif

