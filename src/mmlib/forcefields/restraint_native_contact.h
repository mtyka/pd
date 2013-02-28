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

#ifndef __FFRESTRAINT_NATIVECONTACT_H
#define __FFRESTRAINT_NATIVECONTACT_H

#include "forcefields/forcefield.h" // provides a base class
#include "forcefields/restraintbase.h" // provides a base class

namespace Physics
{


	//-------------------------------------------------
	//
	/// \brief Calculates Native contacts and can apply a constraint to it for umbrella
	/// Type applications etc.
	///
	/// \details 
  ///
  ///    Using this forcefield a peptide or protein can be maintained in its particular backbone 
  ///    conformation using a restraint on the number of native contacts. A form similar to 
	///    Brooks et al. (F. B. Sheinerman and C. L. Brooks III. Calculations on folding of segment b1 of streptococcal protein G. J. Mol. Biol., 278:439456, 1998. ) 
  ///    is used, although the order parameter used here is inverted with respect to theirs (their native state 
  ///    is at $\rho=0.0$ where as in the form used here it is at $\rho=1.0$). The fraction of native 
  ///    contacts in a given conformation is given by
  ///    
  ///    \begin{equation}
  ///      \rho = \frac{1}{N} \sum^N \frac{1}{1 + exp\left[ \vartheta \left( r_{ij} - r_{0} \right) \right]}
  ///    \end{equation}
  ///    
  ///    where the sum runs over all $N$ native contacts, $r_{ij}$ is the inter-atomic distance between the 
  ///    two atoms involved in a given contact while $r_0$ is the distance at which a contact is no longer 
  ///    considered to be formed. $r_0$ was set to 9\AAs. $\vartheta$ is the steepness of the sigmoidal 
  ///    switch function which determines if a given contact is currently formed or not and was set to 
  ///    1.5 giving a comparably smooth transition which 
  ///    appears to work better in a continuous simulation such as Molecular Dynamics.
  ///    Contacts were considered only between C$_\alpha$ atoms and were considered native when their distance 
  ///    was less than $r_{0} - \ln(1.0/0.1 - 1)/\vartheta$ in the starting structure. The extra term  involving 
  ///    the logarithm  ensures that the sigmoidal function is at least 0.9 at the native distance. If this extra 
  ///    term is left out, $\rho$ is significantly smaller than 1.0 at the native structure and applying a restraint 
  ///    to $\rho$ actually distorts the native state since it pushes the longer contacts further together. 
  ///    Note that ``native'' here refers to ``native with respect to a given conformation'', not with 
  ///    respect to the NMR/XRAY structure. 
	///
	/// \author Mike Tyka  
	///
	/// \todo 
	///
	/// \bug 
	///
	class PD_API FF_Restraint_NativeContact : public AtomRestraintForcefieldBase
	{
	public:
		FF_Restraint_NativeContact( WorkSpace &newwspace );

		virtual FF_Restraint_NativeContact* clone() const { return new FF_Restraint_NativeContact(*this); }


		/// prints a little block of parameter information
		virtual void info() const;						

		/// print the list of restraints being used
		virtual void detail() const;				

		/// Cutoff of what's considered a native contact (Angstrom)
		double NativeDist; 

		/// Steepness of sigmoidal switchover of what is considered a contact or not.
		double Steepness; 

		/// p0 is the restraint point, i.e. the proportion of native contacts to be restrained
		/// to. With other words, if p0 is set to 0.5 then the restraint will harmonicaly restrain
		/// the fraction of native contacts to stay around 0.5
		double p0; 

	private:
		virtual void setup();

		int createContacts();

		virtual void calcEnergies();
		virtual void calcForces();

		virtual double calcEnergyAtQ(double newQ);

		/// contact array
		std::vector<IndexPair> contact; 
	};



} // namespace Physics


#endif

