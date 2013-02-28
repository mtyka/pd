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

// Description: Implements various structure pertubators for Monte Carlo Type
// applications, including Normal changes to torsions, RAFT moves etc..

#include "global.h"

#include "basicmoves.h"

#include "maths/rms.h"

#include "library/angleset.h"
#include "library/backbonetorsions.h"
#include "library/rotamerlib.h"

#include "workspace/workspace.h"
#include "workspace/rotbond.h"
#include "workspace/neighbourlist.h"

#include "forcefields/ffparam.h"
#include "forcefields/ffbonded.h"

#include "protocols/minimise.h"
#include "protocols/montecarlo.h"
#include "protocols/dualffminimiser.h"
#include "protocols/md.h"

using namespace Maths;
using namespace Library;
using namespace Physics;
using namespace Protocol;

namespace Manipulator
{	
	int CartesianMove::apply()
	{
		if(fprob>0)
		{
			for(size_t im = 0; im < wspace->atom.size(); im++) 
			{
				if(frand() > fprob)
					continue;
				displaceAtom(im);
			}
		}
		if(ndisplace>0)
		{
			for(size_t i = 0; i < ndisplace; i++)
			{
				size_t im = (size_t)rand() % wspace->atom.size();
				displaceAtom(im);
			}
		}
		return 1;
	}

	void CartesianMove::displaceAtom(size_t i)
	{
		double n1, n2, n3, n4;
		nrand(n1, n2, xsigma);
		nrand(n3, n4, xsigma);
		wspace->cur.atom[i].p.add(n1,n2,n3);
		//printf("move:  %10.7f  %10.7f  %10.7f \n", n1,n2,n3 );
	}


	// -------------------------------------------------------------
	//   Molecule Displacement

	int MoleculeDisplacement::apply(){
		int im, i;

		if(fprob>0){
			for(im = 0; im < wspace->mol.size(); im++) {
				if(frand() > fprob)
					continue;
				displaceMolecule(im);
			}
		}
		if(ndisplace>0){
			for(i=0;i<ndisplace;i++){
				im = rand()%wspace->mol.size();
				displaceMolecule(im);
			}
		}

		return 1;
	};


	void MoleculeDisplacement::displaceMolecule(size_t im){
		double n1, n2, n3, n4, n5, n6, n7, n8;
		nrand(n1, n2, xsigma);
		nrand(n3, n4, xsigma);
		nrand(n5, n6, xsigma);
		nrand(n7, n8, anglesigma);

		dvector disp(n1, n2, n3);
		dvector axis(n4, n5, n6);
		double angle = n7;
		axis.unify();
		matrix3x3 rmat;
		rmat.setToAxisRot(axis, angle);

		int i;
		dvector cog(0, 0, 0);
		for(i = wspace->mol[im].ifirst; i <= wspace->mol[im].ilast; i++) {
			cog.add(wspace->cur.atom[i].p);
		}
		cog.div((double) (wspace->mol[im].ilast - wspace->mol[im].ifirst + 1));
		for(i = wspace->mol[im].ifirst; i <= wspace->mol[im].ilast; i++) {
			wspace->cur.atom[i].p.sub(cog);
			wspace->cur.atom[i].p.mulmat(rmat);
			wspace->cur.atom[i].p.add(cog);
			wspace->cur.atom[i].p.add(disp);
		}
	};

	// -------------------------------------------------------------
	//   Molecule Displacement

	int TIP3P_Move::apply(){
		int im, i;

		if( (movecount % RepeatMoveInterval) == 0 ){
			movedmolecules.clear();
			if(fprob>0){
				for(im = 0; im < wspace->mol.size(); im++) {
					if(frand() > fprob)
						continue;
					applyMoveToMolecule(im);
					movedmolecules.push_back(im);
				}
			}
			if(ndisplace>0){
				for(i=0;i<ndisplace;i++){
					im = rand()%wspace->mol.size();
					applyMoveToMolecule(im);
					movedmolecules.push_back(im);
				}
			}
		}else{
			for(im=0; im<movedmolecules.size(); im++){
				applyMoveToMolecule(im);
			}
		}
		movecount++;
		return 1;
	};


	void TIP3P_Move::applyMoveToMolecule(size_t im){
		double n0,n1, n2, n3, n4, n5, n6, n7, n8,n9;
		nrand(n1, n2);
		nrand(n3, n4);
		nrand(n5, n6);
		nrand(n7, n8);
		nrand(n9, n0);

		dvector disp(n1, n2, n3);
		disp.mul(trans_sigma);
		dvector axis(n4, n5, n6);
		dvector noise;
		axis.unify();
		matrix3x3 rmat;
		rmat.setToAxisRot(axis, n7 * rot_sigma);

		int i;
		dvector cog(0, 0, 0);
		for(i = wspace->mol[im].ifirst; i <= wspace->mol[im].ilast; i++) {
			cog.add(wspace->cur.atom[i].p);
			wspace->atom[i].setMoved(true);        // flag that this atom has been changed (it will change shortly)
		}
		cog.div((double) (wspace->mol[im].ilast - wspace->mol[im].ifirst + 1));
		for(i = wspace->mol[im].ifirst; i <= wspace->mol[im].ilast; i++) {
			wspace->cur.atom[i].p.sub(cog);
			wspace->cur.atom[i].p.mulmat(rmat);
			wspace->cur.atom[i].p.add(cog);
			wspace->cur.atom[i].p.add(disp);
		}

		int first = wspace->mol[im].ifirst;
		wspace->cur.atom[first+0].p.add(n1*bond_sigma,n6*bond_sigma,n7*bond_sigma);
		wspace->cur.atom[first+1].p.add(n2*bond_sigma,n5*bond_sigma,n8*bond_sigma);
		wspace->cur.atom[first+2].p.add(n3*bond_sigma,n4*bond_sigma,n9*bond_sigma);

		/*
		int first = wspace->mol[im].ifirst;
		dvector bond;
		bond.diff(wspace->cur.atom[first+1].p,wspace->cur.atom[first].p);
		bond.unify();
		bond.mul(n9);
		wspace->cur.atom[first+1].p.add(bond);
		bond.diff(wspace->cur.atom[first+2].p,wspace->cur.atom[first].p);
		bond.unify();
		bond.mul(n0);
		wspace->cur.atom[first+2].p.add(bond);
		*/


	};

	// -------------------------------------------------------------
	//   Molecule Displacement

	// changes each rotatable torsion by a normal function
	int NormalTorsionalMove::apply(){
		int i;
		double nrandnum1;
		double nrandnum2;
		int dihedrand;
		double change;
		double changeproportion = 0.15;
		double frq_backbonejump = 0.005;

		int moveseverity = 0;

		RotBond& rotbond = wspace->rotbond();

		for(i = 0; i < rotbond.size(); i++) {
			if( (rotbond[i].Type == RotatableBond::Double) ||
				(rotbond[i].Type == RotatableBond::Omega) ||
				(rotbond[i].Type == RotatableBond::Planar)) {

					// do nothing for these types....

			} else if((rotbond[i].Type == RotatableBond::Phi) || (rotbond[i].Type == RotatableBond::Psi)) {

				if(frand() < frq_backbonejump) {
					nrand(nrandnum1, nrandnum2, DegToRad(251.0) ); //120
					rotbond.rotate(i, nrandnum1);
					moveseverity = max(moveseverity, 4);
				} else {
					nrand(nrandnum1, nrandnum2, DegToRad(3.90) ); //15
					rotbond.rotate(i, nrandnum1);
					moveseverity = max(moveseverity, 1);
				}

			} else {

				// is is a non phi or psi single bond ....

				if(frand() < changeproportion) {
					nrand(nrandnum1, nrandnum2, DegToRad(6.9) ); //20
					dihedrand = rand() % 3 - 1;
					change = (double)dihedrand * DegToRad(120.0) + nrandnum1;
					rotbond.rotate(i, change);
				} else {
					nrand(nrandnum1, nrandnum2, DegToRad(3.9) ); //15
					rotbond.rotate(i, nrandnum1);
				}
				moveseverity = max(moveseverity, 1);
			}
		}
		return moveseverity;
	};




	// changes each rotatable torsion by a normal function
	int SidechainTorsionalMove::apply(){
		int i;
		double nrandnum1;
		double nrandnum2;
		int dihedrand;

		int moveseverity = 0;

		RotBond& rotbond = wspace->rotbond();

		for(i = 0; i < rotbond.size(); i++) {
			if((rotbond[i].Type == RotatableBond::Single)) {
				if(frand() < p120) {
					dihedrand = rand() % 3 - 1;
					rotbond.rotate(i, dihedrand * DegToRad(120.0) );
				}
				if(frand() < pnormal) {
					nrand(nrandnum1, nrandnum2, sdnormal); //15
					rotbond.rotate(i, nrandnum1);
				}
				moveseverity = max(moveseverity, 1);
			}
		}
		return moveseverity;
	};


	// changes each backbone (psi/psi) torsion by a normal function
	int BackboneTorsionalMove::apply(){
		int i;
		double nrandnum1;
		double nrandnum2;
		int moveseverity = 0;

		RotBond& rotbond = wspace->rotbond();

		for(i = 0; i < rotbond.size(); i++) {
			if((rotbond[i].Type == RotatableBond::Phi) || (rotbond[i].Type == RotatableBond::Psi)) {
				if(frand() < pnormal) {
					moveseverity = 1; // removed max(moveseverity, 1); : because will always be 1 if we have any cases of (frand() < pnormal)
					nrand(nrandnum1, nrandnum2, sdnormal); //15
					if(maxresidues > 0)
						rotbond.rotate(i, nrandnum1, maxresidues);
					else
						rotbond.rotate(i, nrandnum1);
				}
			}
		}
		return moveseverity;
	};


	int BackbonePropagationMove::apply(){
		int ir;
		int relative;
		double phi, psi;
		int changes = 0;

		for(ir = 1; ir < (wspace->res.size() - 1); ir++) {
			if(frand() > pprop)
				continue;

			wspace->calcResiduePhiPsi(ir, phi, psi);
			relative = (rand() % 2) * 2 - 1; // +1 / -1
			wspace->setResiduePhiPsi(ir + relative, psi, psi);

			changes++;
		}
		if(changes == 0)
			return 0;
		return 3;
	};


	int CartesianBlockMove::apply(){
		if(frand() > fglobal)
			return 0;

		size_t startir;
		size_t endir;
		int i;
		double n1, n2, n3, n4, n5, n6, n7, n8;


		startir = rand() % wspace->res.size();
		endir = min((int)(wspace->res.size() - 1), int(startir + rand() % 2 + 1));

		if(endir > (wspace->res.size() - 1))
			endir = (wspace->res.size() - 1);
		nrand(n1, n2, xsigma);
		nrand(n3, n4, xsigma);
		nrand(n5, n6, xsigma);
		nrand(n7, n8, anglesigma);

		dvector disp(n1, n2, n3);
		dvector axis(n4, n5, n6);
		double angle = n7;
		axis.unify();
		matrix3x3 rmat;
		rmat.setToAxisRot(axis, angle);

		dvector cog(0, 0, 0);
		for(i = wspace->res[startir].ifirst; i < wspace->res[endir].ilast; i++) {
			cog.add(wspace->cur.atom[i].p);
		}
		cog.div((double) (i - wspace->res[startir].ifirst));
		for(i = wspace->res[startir].ifirst; i < wspace->res[endir].ilast; i++) {
			wspace->cur.atom[i].p.sub(cog);
			wspace->cur.atom[i].p.mulmat(rmat);
			wspace->cur.atom[i].p.add(cog);
			wspace->cur.atom[i].p.add(disp);
		}
		return 1;
	};



	int PeptideGroupMove::apply(){
		if(frand() > fglobal)
			return 0;

		int ir;
		dvector axis;
		matrix3x3 rmat;
		double angle, n2;
		int iCA0, iC0, iO0, iN1, iH1, iCA1;

		for(ir = 0; ir < (wspace->res.size() - 1); ir++) {
			if(frand() > prop)
				continue;

			if((wspace->res[ir + 1].letter == 'p') || (wspace->res[ir + 1].letter == 'P'))
				continue;
			iCA0 = wspace->res[ir].iCA;
			iCA1 = wspace->res[ir + 1].iCA;

			if((iCA0 < 0) || (iCA1 < 0))
				continue; // at least the CAs are required

			iC0 = wspace->res[ir].iC;
			iO0 = wspace->res[ir].iO;
			iN1 = wspace->res[ir + 1].iN;
			iH1 = wspace->res[ir + 1].iH;

			nrand(angle, n2, anglesigma);
			axis.diff(wspace->cur.atom[iCA0].p, wspace->cur.atom[iCA1].p);
			rmat.setToAxisRot(axis, angle);

			if(iC0 >= 0) {
				wspace->cur.atom[iC0].p.sub(wspace->cur.atom[iCA0].p);
				wspace->cur.atom[iC0].p.mulmat(rmat);
				wspace->cur.atom[iC0].p.add(wspace->cur.atom[iCA0].p);
			}
			if(iO0 >= 0) {
				wspace->cur.atom[iO0].p.sub(wspace->cur.atom[iCA0].p);
				wspace->cur.atom[iO0].p.mulmat(rmat);
				wspace->cur.atom[iO0].p.add(wspace->cur.atom[iCA0].p);
			}
			if(iN1 >= 0) {
				wspace->cur.atom[iN1].p.sub(wspace->cur.atom[iCA0].p);
				wspace->cur.atom[iN1].p.mulmat(rmat);
				wspace->cur.atom[iN1].p.add(wspace->cur.atom[iCA0].p);
			}
			if(iH1 >= 0) {
				wspace->cur.atom[iH1].p.sub(wspace->cur.atom[iCA0].p);
				wspace->cur.atom[iH1].p.mulmat(rmat);
				wspace->cur.atom[iH1].p.add(wspace->cur.atom[iCA0].p);
			}


		}

		return 1;
	};





	int BackbonePhiPsiSetMove::apply(){
		int ir;
		int changes = 0;

		dvector *oldpos = new dvector[wspace->atom.size()];
		for(int ii = 0; ii < wspace->atom.size(); ii++)
			oldpos[ii].setTo(wspace->cur.atom[ii].p);

		int changable[100];
		int nchangeable = 0;
		int c;

		if(frand() < (fchange)) {
			ir = rand() % (wspace->res.size() - simultaneous + 1);

			changeConsResiduesRandomly(ir, simultaneous);
			changes++;
			for(c = 0; c < wspace->res.size(); c++) {
				if(c == ir)
					c += simultaneous;
				changable[nchangeable] = c;
				nchangeable++;
			}

			if(local != 0) {
				if((ir > endlength) && ((ir + simultaneous) < (wspace->res.size() - endlength)));
				else
					nchangeable = 0;

				nchangeable = 0;
			}

			if(sidechainfix != 0) {
				// sort out any sterically challanged sidechains
				SidechainRotamerLibMove scrlp(*wspace, *rotlib, 0.0, 6.0);
				scrlp.apply();
			}
		}
		// printf("Changes: %d \n", changes);
		delete[]oldpos;

		return changes;
	}


	void BackbonePhiPsiSetMove::changeResidueRandomly(int ir, bool weighted)
	{
		const BackboneTorsionLibrary& set = angset->getBackboneTorsionSet( wspace->res[ir].letter );

		conformer_type angcount = (conformer_type)set.size();	
		conformer_type anglepair;
		if( weighted )
		{
			// weighted probability
			double *problist = new double[angcount];
			for(conformer_type p = 0; p < angcount; p++)
				problist[p] = set.getPropensity(p);
			anglepair = chooseByProbability(problist, angcount);
			delete[] problist;
		}
		else
		{
			// uniform probability
			anglepair = (conformer_type)(rand()%angcount);
		}

		wspace->setResiduePhiPsi(ir, set.getPhi(anglepair), set.getPsi(anglepair));	
	}

	void BackbonePhiPsiSetMove::changeConsResiduesRandomly(int sir, int nir){
		// changes consecutive residues from sir to sir+nir
		int ir;
		for(ir = sir; ir < (sir + nir); ir++) {
			if(ir < wspace->res.size())
				changeResidueRandomly(ir);
		}
	}


	BlockRotationMove::BlockRotationMove(
		WorkSpace& _wspace,
		Physics::Forcefield& _stericff,
		double _fchange)
		: MoveBase(_wspace) 
	{
		stericff = &_stericff;
		fchange = _fchange;
	}

	BlockRotationMove::~BlockRotationMove()
	{
	}

	int BlockRotationMove::apply(){
		using namespace Protocol;

		if(frand() > (fchange))
			return 0;

		int i;
		int Steps = 6;
		int Step = 0;
		int minsteps = 100;
		double stepangle = DegToRad( 10.0 );
		SnapShot *station = new SnapShot[Steps * 2 + 1];
		int nstations = 0;
		int seglength = rand() % 5 + 4;
		int ir1 = rand() % (wspace->res.size() - seglength);
		int ir2 = ir1 + seglength;

		if(ir2 >= (wspace->res.size() - 2))
			ir2 = wspace->res.size() - 1;
		if(ir1 <= 1)
			ir1 = 0;

		printf("Rotating %d -> %d \n", ir1, ir2);

		DualFFMinimiser minimise(*stericff);
		minimise.Steps = minsteps;
		minimise.SlopeCutoff = 0.05 * PhysicsConst::kcal2J / PhysicsConst::Na;

		station[nstations] = wspace->save();
		nstations++;

		for(Step = 0; Step < Steps; Step++) 
		{
			rotateBlock(ir1, ir2, stepangle);
			minimise.runcore();
			station[nstations]= wspace->save();
			nstations++;
		}

		// restore to starting position
		wspace->load(station[0]);

		for(Step = 0; Step < Steps; Step++) {
			rotateBlock(0, 8, -stepangle);
			minimise.runcore();

			station[nstations] = wspace->save();
			nstations++;
		}

		int lowestStation;
		double lowestStationEnergy = station[1].epot;
		for(i = 1; i < nstations; i++) {
			if(double (station[i].epot) < lowestStationEnergy) {
				lowestStationEnergy = station[i].epot;
				lowestStation = i;
			}
		}

		wspace->load(station[lowestStation]);
		delete[]station;


		return 1;
	}




	int BlockRotationMove::rotateBlock(int ir1, int ir2, double angle){
		int c;


		dvector axis;
		axis.diff(wspace->cur.atom[wspace->res[ir1].iCA].p, 
			wspace->cur.atom[wspace->res[ir2].iCA].p);
		dvector origin;
		origin.setTo(wspace->cur.atom[wspace->res[ir1].iCA].p);

		for(c = wspace->res[ir1].iCA; c < wspace->res[ir2].iCA; c++) {
			wspace->cur.atom[c].p.rotateAxis(angle, origin, axis);
		}

		return 0;
	}




	SidechainRotamerLibMove::SidechainRotamerLibMove(
		WorkSpace& _wspace,
		const Library::RotamerLibrary& _rotlib,
		double _pprop,
		double _forcedchangecutoff, 
		RotamerMode _mode)
		: RotamerApplicatorBase(_wspace,_rotlib, _mode) 
	{
		pprop = _pprop;
		forcedchangecutoff = _forcedchangecutoff;
	}

	int SidechainRotamerLibMove::apply(){
		int ir, irot;
		int nrots;
		int nvalidrots;
		double rotamersterics[100];
		int validrotamer[100];
		double stericlimit = 4.0;
		bool changerotamer;

		for(ir = 0; ir < wspace->res.size(); ir++) {
			changerotamer = false;
			if(frand() < pprop)
				changerotamer = true;

			nrots = getRotSterics(ir, &rotamersterics[0], changerotamer, forcedchangecutoff);

			if((!changerotamer) && (nrots <= 1))
				continue; // dont do anything if rotamer is ok and we dont force a change

			if(nrots <= 1)
				continue;// no rotamers to chose from
			// if(!changerotamer) printf("Forced change of rotamer ir=%d : steric score: %lf -> choices %d\n",
			// ir,rotamersterics[0],nrots);
			nvalidrots = 0;
			for(irot = 1; irot < nrots; irot++) {
				// printf("%lf \n",rotamersterics[irot]);
				if(rotamersterics[irot] < stericlimit) {
					validrotamer[nvalidrots] = irot;
					nvalidrots++;
				}
			}
			if(nvalidrots <= 0) { // no rotamers to chose from, at least chose the one with the smallest steric score
				double lowestSteric = DBL_MAX;
				int lowestStericIndex = 0;

				for(irot = 0; irot < nrots; irot++) {
					// printf("%lf \n",rotamersterics[irot]);
					if(rotamersterics[irot] < lowestSteric) {
						lowestSteric = rotamersterics[irot];
						lowestStericIndex = irot;
					}
				}
				if(lowestStericIndex == 0)
					continue; // dont do anything if the current conformation is already the lowest steric
				irot = lowestStericIndex; // otherwise change sidechain to the rotamer witht he lowest steric value
			} else { // if several rotamers are fine choose randomly amongst the ok ones
				irot = validrotamer[rand() % nvalidrots];
			}

			// printf("Changing: %d \n",irot);
			//rotlib->changeSidechainConformation(wspace, ir, irot - 1);
			applyRotamer(ir,irot-1);
			// wspace->outtra.append();
		}


		return 1;
	}


	// This function will do two related things:
	// a) it will check the steric situation of the current sidechain conformation or the
	// residue ir and suggest a change if above a certain threshold (scorecutoff) (using b) )
	// b) if tryallrotamers = true it will assess the entire rotamer distribution

	int SidechainRotamerLibMove::getRotSterics(int ir, double *rotsterics, bool tryallrotamers, double scorecutoff){
		dvector *saveposition = NULL;
		int irot = 0;
		int i, j;
		double sqrdistij;
		double radius = 3.4;
		double sqrradius = sqr(radius);
		double score;
		double scorei;
		double scoreij;

		// first save atom positions
		saveposition = new dvector[wspace->res[ir].ilast - wspace->res[ir].ifirst + 1];
		for(i = wspace->res[ir].ifirst; i <= wspace->res[ir].ilast; i++) 
		{
			saveposition[i - wspace->res[ir].ifirst].setTo(wspace->cur.atom[i].p);
		}

		// got through every possible rotamer, apply it, and calculate the
		// local steric impact
		// the first round through this loop accesses the current rotamer, its steric situation
		// information will be saved in rotsterics[0]
		
		// } while(rotlib->changeSidechainConformation(wspace, ir, irot - 1) == 0);
		for( int k = 0; k < nRot(ir)+1; k++ )
		{
			if( k != 0 ) applyRotamer(ir,k-1);

			// estimate the atomic overlaps;
			score = 0;
			for(i = wspace->res[ir].ifirst; i <= wspace->res[ir].ilast; i++) {
				if(wspace->atom[i].Z == 1)
					continue; // ignore hydrogens
				if(i == wspace->res[ir].iCA)
					continue; // skip backbone atoms
				if(i == wspace->res[ir].iC)
					continue; // skip backbone atoms
				if(i == wspace->res[ir].iN)
					continue; // skip backbone atoms
				if(i == wspace->res[ir].iH)
					continue; // skip backbone atoms
				if(i == wspace->res[ir].iO)
					continue; // skip backbone atoms
				scorei = 0;
				for(j = 0; j < wspace->atom.size(); j++) {
					// skip self interactions
					if((j >= wspace->res[ir].ifirst) && (j <= wspace->res[ir].ilast)) {
						j = wspace->res[ir].ilast + 1;
						continue;
					}
					if(wspace->atom[j].Z == 1)
						continue; // ignore hydrogens

					// now calculate the squaredistance between i & j

					sqrdistij = SqrDistanceBetween(
						wspace->cur.atom[i].p.x,
						wspace->cur.atom[i].p.y,
						wspace->cur.atom[i].p.z,
						wspace->cur.atom[j].p.x, 
						wspace->cur.atom[j].p.y, 
						wspace->cur.atom[j].p.z);
					if(sqrdistij > sqrradius)
						continue;
					scoreij = sqr(1.0 - (sqrdistij / sqrradius));
					scorei += scoreij;
				}
				// printf("scorei: %lf \n",scorei);
				score += scorei;
			}
			//printf("Rot: %3d score = %lf \n",irot,score);
			//wspace->ene.epot = score * PhysicsConst::kcal2J / PhysicsConst::Na;
			//if(score < 3.5) wspace->outtra.append();

			rotsterics[irot] = score;
			irot++;
			if((irot == 1) && (!tryallrotamers) && // if the current conformation (irot==1) is ok sterically
				(score < scorecutoff)) { // and caller has not explicitly asked for all rotameters, we quit
					break;
			}
			// irot = 0 means current conformation, irot = 1 means the 0th rotamer,
			// hence the -1 in the line below
		}

		// restore positions
		for(i = wspace->res[ir].ifirst; i <= wspace->res[ir].ilast; i++) 
		{
			wspace->cur.atom[i].p.setTo(saveposition[i - wspace->res[ir].ifirst]);
		}
		delete[]saveposition;

		return irot;
	}
}


