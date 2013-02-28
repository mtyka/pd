PD – A Molecular Mechanics Engine for Python

PD is a free, modular C++ library for biomolecular simulation with a 
flexible and scriptable Python interface. 

Copyright (C) 2003-2013 Mike Tyka and Jon Rea

Contact: mike.tyka@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.


PD is a highly modular, cross-platform molecular mechanics suite for Python. 
The library design is focused on providing unprecedented flexibility, modularity 
and adaptability without compromising speed of execution. The core algorithms 
are written in C++ for maximum efficiency and make extensive use of 
object-oriented design patterns. All features are automatically wrapped into 
Python using SWIG, a widely used wrapper generator, circumventing the need for 
the developer to know anything about the exporting process while still allowing 
Python to be used as the input script language. This allows the user to set up 
simple simulations quickly and easily while providing powerful language features 
for the more advanced user. All abstract concepts in molecular simulations such 
as force fields and simulation methods are represented by classes which interact 
with one another in a highly modular fashion, allowing the user to create new 
algorithms without coding a single line of C++. For the developer the object 
oriented design of the core library allows easy reuse of code and rapid development 
and testing of new algorithms hence accelerating scientific progress. Furthermore, 
it is possible to write separate modules which interact with the core library but 
are separately compiled and maintained and enable scientific groups to provide 
new algorithms to the code base without the need for code patches to the core 
library. At this early stage, the library already provides a large number of 
widely used features: Various molecular dynamics protocols (e.g. Replica Exchange 
Simulation) and Monte Carlo algorithms as well as methods for free-energy 
calculation, quasi-harmonic and normal mode analysis, have been implemented; 
Algorithms for homology modeling and docking are under development. Commonly 
used force fields AMBER, CHARMM and OPLS are currently supported, including 
periodic boundary simulations. Implicit solvation, such as surface area based 
models and Generalized-Born solvent models are also implemented.


FEATURES


Written in C++, wrapped into Python

Exposure into Python allows very simple as well as highly complex use of the 
features giving the user a maximum of flexibility when setting up simulations.
The design allows classes to recombine into novel algorithms at runtime (i.e. from python).
Interface generated automatically using SWIG: No need to write interface code 
after a new feature is added.

Forcefield functional forms: 

Bonded: Harmonic Bonds & Angles, Fourier & Harmonic Torsions & Impropers
NonBonded: Lennard Jones Van der Waals, Simple Electrostatics. Energy and Force switching.
Spaces: Periodic boundary conditions (Orthogonal)
Solvation: Distance-dependent Dielectric Constant, Solvent Accessible Surface 
           Area (SASA) based solvation, Generalised Born / Surface Area (GB/SA)
Knowledge Based: Residue-Residue Interaction matrix based
Restraints: Harmonic Cartesian, Distance & Torsion Restraints & various other custom forces.
Other: Go forcefield, Softcore VdW Forcefield
OpenMP allows parallelism of NonBonded calculations across the CPUs of Multicore/MultiCPU machines.

Forcefields: 

AMBER: (ff94, ff96, ff99, ff03),
CHARMM: (Charmm19, Charmm22),
Other: OPLS/aa, Mizyawa & Jernigan, RAFT

Protocols:

Gradient based Minimisation (Steepest Descent, Conugate Gradients & Torsional)
Monte Carlo, Monte Carlo with Minimiation including a large variety of Structure 
Perturbation (Move Types) including Rotamer Libraries & Fragment Insertions
Conformational Space Annealing (CSA)
Homology Modelling: Threading & Loop Builder
Molecular Dynamics: Beeman , Verlet & Velocity Verlet Integrators, Andersen and 
Berendsen Thermostats
Langevin Dynamics
Replica Exchange Dynamics (REMD/REMC)
Numerical Second Derivate Matrix (Hessian) Calculation and Normal Mode Analysis
Quasi Harmonic Analysis
Free Energy Calculations: Free Energy Pertubation , Thermodynamic Integration, 
Umbrella Sampling

File Formats

Input: PDB, TRA
Output: PDB, DCD/PSF, TRA



INSTALLATION 
See doc/installation/install.html for detailed information on installation
and/or compilation of the code.

DOCUMENTATION
See doc/index.html for more detailed documenation and tutorials


DIRECTORY STRUCTURE

doc/   Documentation
src/   Source Code and GNU Makefile
bin/   Binaries (for source distributions only after compilation)
lib/   Forcefield and other run-time libraries
vc/    Visual C (MS VisualStudio 2005) project files for Windows


Disclaimer:

* THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

