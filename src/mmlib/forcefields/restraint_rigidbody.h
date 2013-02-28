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


#ifndef __FF_RESTRAINT_RB_H
#define __FF_RESTRAINT_RB_H

#include "system/fundamentals.h"
#include "forcefields/forcefield.h" // Provides the base class
#include "forcefields/restraintbase.h" // provides a base class

#include "pickers/pickbase.h"          // provides a member variable
#include "workspace/pospointer.h"   // provides a member variable
#include "workspace/workspace.fwd.h"

class PD_API Particle;

namespace Physics
{


//-------------------------------------------------
//
/// \brief  Implements a special Rigid Body Restraint which does not affect the internal motion of
///          the set of atoms involved. This is merely the base class for three sepearte classes.
///
/// \details 
///    In free energy calculations and other situations it can be advantageous to be able to restrain the rigid body
///    motion of a set of atoms in such a way that the free energy of the restraint is analytically computable.
///    In J. Hermans and L.Wang JACS 1997, 119, 2707-2714 a set of such restraints is described.
///    
///    Three restraints are applied; A positional restraint preventing movement from the centre of the simulation 
///    box (3 degrees of freedom), an angular restraint maintaining an orientation along a given vector (2 degrees 
///    of freedom) and a dihedral restraint preventing rotation around that vector (1 degree of freedom). The free 
///    energies associated with these 3 restraints are $\Delta A_x$,$\Delta A_\theta$ and $\Delta A_\chi$. For the 
///    sake of completeness it ought to be mentioned that a further free- energy component ($\Delta A_s$) would be 
///    present for symmetric molecules.
///    
///    More information can be found in J. Hermans and L.Wang JACS 1997, 119, 2707-2714i
///
///    The three restraints are implemented as seperate child classes of this abstract base class and are
///    intended to be used in conjunction.
///
/// \author Mike Tyka 
///
/// \todo  
///
/// \bug 
///
	class PD_API FF_BodyRestraint: public AtomRestraintForcefieldBase{
	public:

		FF_BodyRestraint(WorkSpace &newwspace): 
			AtomRestraintForcefieldBase(newwspace) 
		{
			k = 0;
			Power = 2;
			MassWeight = false;
		}
		virtual FF_BodyRestraint* clone() const = 0;

		double Power; ///< Power of the restraint (epot = k*|x-x0|^Power )
		
		bool MassWeight;
		

	protected:

		virtual void setup(){
			RestraintForcefieldBase::setup();
		}

	};







//-------------------------------------------------
//
/// \brief  Translational Rigid Body Restraint - restrains the center of mass of a arbitrary set of atoms
///
/// \details 
///    The position restraint acts on the centre of mass of the solute. The potential energy
///    is 
///    
///    \begin{equation}
///      \mathscr{V}(\mathbf{x}) = \frac{1}{2} k_\mathbf{x} \left(\mathbf{x} - \mathbf{x}_0 \right)^2
///    \end{equation}
///    where $k_\mathbf{x}$ is the translational restraint constant.
///    
///    The free energy of applying the restraint is 
///    \begin{eqnarray}
///    \Delta A_x & = & -k_BT \ln \left[1/V_0 \int_{V_0} exp(- \mathscr{V}(\mathbf{x})/k_BT) \mathrm{d} \mathbf{x} \right] \\
///     & \approx & k_BT \ln V_0 - k_BT \ln \left( \frac{2\pi k_BT}{k_\mathbf{x}} \right)^{2/3}  
///    \end{eqnarray}
///    where $V_0$ is the volume of the simulation box, $T$ is the temperature and $k_B$ is the Boltzmann constant. 
///    The last approximation is valid when the size of the simulation box is much larger than the extent of the space 
///    accessible to the restrained molecule, or more formally when $d\sqrt{\frac{k_\mathbf{x}}{k_BT} }$ is greater 
///    than $ \approx 10$, where d is the shortest length of the simulation box. We used a value of $k_\mathbf{x} = 50$ \kunit, 
///    thus this condition is fulfilled  by a wide margin.
///    
/// \author Mike Tyka 
///
/// \todo 
///
/// \bug 
///
	class PD_API FF_BodyRestraint_Position:public FF_BodyRestraint{
	public:

		FF_BodyRestraint_Position( WorkSpace &newwspace ): FF_BodyRestraint(newwspace) {
			name = "BodyRestraint(Pos)";
			ShortName = "RestPos";
			k = 0;
			Power = 2;
		}
		virtual FF_BodyRestraint_Position* clone() const { return new FF_BodyRestraint_Position(*this); }
	
		double Power; ///< Power of the restraint (epot = k*|x-x0|^Power )

		virtual double calcEnergyAtQ(double newQ);
		virtual void info() const; ///< prints a little block of parameter information  

	protected:
		virtual void setup();
		virtual void  calcForces();
	};

	enum FunctionalForm { Harmonic, Cosine };






//-------------------------------------------------
//
/// \brief  Angular Rigid Body Restraint - restrains the angule between the centers of mass of two arbitrary set of of atoms and a chosen absolute vector 
///
/// \details
///  Angular Rigid Body Restraint - restrains the angule between the centers of mass of two arbitrary set of of atoms and a chosen absolute vector 
///
///    The angular restraint acts on the angle between the vector $\mathbf{r}$, connecting the centre of masses $\mathbf{x}_i$ and 
///    $\mathbf{x}_j$ of two sets of atoms within the molecule, and a reference vector $\mathbf{e}_\theta$.
///    
///    \begin{equation}
///      \mathscr{V}(\mathbf{\theta}) = \frac{1}{2} k_\mathbf{\theta} \left(1 - cos \theta \right)
///    \end{equation}
///    where $\theta$ is the angle between the two vectors and $k_\mathbf{\theta}$ is the associated force constant.
///    
///    The free energy of applying this angular restraint is 
///    \begin{equation}
///    \Delta A_\theta = -k_B T \ln \left\{ (k_BT/k_\theta)  
///    \left[ 1 - exp(-k_\theta / k_BT ) \right] \right\}
///    \end{equation}
///    
///
/// \note
///    For the angular and dihedral restraints we have implemented harmonic as well as
///    cosine terms;
///
///
/// \author Mike Tyka 
///
/// \todo 
///
/// \bug 
///
	class PD_API FF_BodyRestraint_Angle:public FF_BodyRestraint{
	public:

		FF_BodyRestraint_Angle( WorkSpace &newwspace ): FF_BodyRestraint(newwspace) {
			name = "BodyRestraint(Angle)";
			ShortName = "RestAng";
			m_AngleTerm_Vector.setTo(0,0,1);
			Form = Harmonic;
		}

		virtual FF_BodyRestraint_Angle* clone() const { return new FF_BodyRestraint_Angle(*this); }
		virtual double calcEnergyAtQ(double newQ);
		virtual void info() const; ///< prints a little block of parameter information  
		FunctionalForm Form;
	protected:
		Maths::dvector m_AngleTerm_Vector;
		virtual void  setup();
		virtual void  calcForces();
	};







//-------------------------------------------------
//
/// \brief  Torsional Rigid Body Restraint - restrains the normal vector of the plane 
///         connecting three center of masses of three arbitrary groups of atoms with a reference vector
///
/// \details 
///    The final degree of freedom is the dihedral angle $\chi$ defined by two planes. The first plane goes 
///    through the points $\mathbf{x}_i$, $\mathbf{x}_j$ and a third point $\mathbf{x}_k$.
///    The second plane goes through $\mathbf{x}_i$, $\mathbf{x}_j$ and a reference vector $\mathbf{e}_\chi$.
///    
///    
///    \begin{equation}
///      \mathscr{V}(\mathbf{\chi}) = \frac{1}{2} k_\mathbf{\chi}  \chi^2 
///    \end{equation}
///    
///    The free energy of the restraint is 
///    
///    \begin{equation}
///      \Delta A_\chi \approx - \frac{1}{2} k_BT \ln \left( 2\pi k_BT / k_\mathbf{\chi} \right) \: \: \: \:\mathrm{ for } \: \: \: \: k_\chi \gg k_BT
///    \end{equation}
///    
///
/// \author Mike Tyka 
///
/// \todo 
///
/// \bug 
///
	class PD_API FF_BodyRestraint_Dihedral:public FF_BodyRestraint{
	public:

		FF_BodyRestraint_Dihedral( WorkSpace &newwspace ): FF_BodyRestraint( newwspace ) {
			name = "BodyRestraint(Dihedral)";
			ShortName = "RestDihe";
			m_DihedralTerm_Vector.setTo(1,0,0);
			Form = Harmonic;
		}
		virtual FF_BodyRestraint_Dihedral* clone() const { return new FF_BodyRestraint_Dihedral(*this); }

		virtual double calcEnergyAtQ(double newQ);
		virtual void info() const; ///< prints a little block of parameter information  
		FunctionalForm Form;

	protected:
		Maths::dvector m_DihedralTerm_Vector;
		virtual void  setup();
		virtual void  calcForces();
	};

} // namespace Physics

#endif

