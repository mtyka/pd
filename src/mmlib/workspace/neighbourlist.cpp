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

#include "space.h"
#include "workspace.h"
#include "bondorder.h"

#include "neighbourlist.h"

using namespace Maths;

NeighbourListBase::NeighbourListBase() :
WorkSpaceComponentBase(NULL),
fnbor(NULL),
neighborlistspace(NULL),
currentMaxNeighbors(0),
Enabled(true),
Cutoff(0.0), 
Padding(0.5),
CalcShadow(false)
{
	m_FullUpdateCount = 0; 
}

NeighbourListBase::NeighbourListBase(const NeighbourListBase &copy)
{
	(*this) = copy;
}

NeighbourListBase::~NeighbourListBase()
{
	clear();
}

void NeighbourListBase::clear()
{
	delete[] fnbor;
	delete[] neighborlistspace;
}

void NeighbourListBase::reinit( WorkSpace* _wspace )
{
	wspace = _wspace;
	reassignList();
}


void NeighbourListBase::requestCutoff(double newcutoff)
{
	if(newcutoff > Cutoff) Cutoff = (double)newcutoff;
	wspace->boundary().setupCellDimensions(Cutoff);
}

void NeighbourListBase::setPadding( double newpadding )
{
	if(newpadding < 0) printf("WARNING Neighborlist Padding cannot be < 0\n");
	Padding = newpadding;
}

void NeighbourListBase::reassignList()
{
	printf("Allocating new full neighbor list\n");
	int nlist_natom = wspace->atom.size();

	delete[] fnbor;
	fnbor = new NeighbourData[nlist_natom];

	if(fnbor == NULL) {
		throw(OutOfMemoryException("Neighbour list allocation failed - out of memory!"));
	}

	reserveMemoryFor(Maths::minmax(4,nlist_natom,100));
}


size_t NeighbourListBase::memuse(int level)
{

	unsigned neighborlistsize = currentMaxNeighbors * (sizeof(int) + sizeof(double) + sizeof(int));
	unsigned bytesrequired = sizeof(*this) + wspace->atom.size() * neighborlistsize;
	if(level>1)printf("    ptr_nlist: %d(%d) --> %3.2lf Mb \n",
		wspace->atom.size(),currentMaxNeighbors,
		double(bytesrequired)/(1024.0*1024.0));
	return bytesrequired;
}

void NeighbourListBase::reserveMemoryFor(int newmaxneighbors)
{

	//if( newmaxneighbors > wspace->atom.size() ) newmaxneighbors = wspace->atom.size(); // There is never a need for '>nlist_natom'
	if( newmaxneighbors < currentMaxNeighbors ) return; // If we already have enough memory, there is no need to reallocate...

	// Free the current memory allocation
	delete[] neighborlistspace;

	unsigned bytesrequired;
	unsigned neighborlistsize;

	neighborlistsize = newmaxneighbors * (sizeof(int) + sizeof(int));
	bytesrequired = wspace->atom.size() * neighborlistsize;
	printf("Growing Neighborlist: Atoms = %d  MaxNeighbors/Atom = %d  Memory use = %3.2lf Mb \n",
		wspace->atom.size(),newmaxneighbors,double(bytesrequired)/(1024.0*1024.0));

	neighborlistspace = new char[bytesrequired];
	for(int i = 0; i < wspace->atom.size(); i++) {
		fnbor[i].n = 0;
		fnbor[i].i = (int *) (&neighborlistspace[i * neighborlistsize]);
		fnbor[i].Type = (int *)
			(&neighborlistspace[i * neighborlistsize + newmaxneighbors * (sizeof(int))]);
	}

	currentMaxNeighbors = newmaxneighbors;


	// clear neighbor arrays
	for(size_t i=0;i<wspace->atom.size();i++) fnbor[i].n=0;

	return;
}

void NeighbourList_32Bit_Base::reassignList()
{
	printf("Allocating new full neighbor list\n");
	int nlist_natom = wspace->atom.size();

	delete[] fnbor;
	fnbor = new NeighbourData[nlist_natom];

	if(fnbor == NULL) {
		throw(OutOfMemoryException("Neighbour list allocation failed - out of memory!"));
	}

	reserveMemoryFor(nlist_natom*Maths::minmax(4,nlist_natom,100));
}

// in this incarnation of the neighbor list newmaxneighbors and maxneighbors
// mean the total number of neighbors - NOT the number of neihbors per atom! 
void NeighbourList_32Bit_Base::reserveMemoryFor(int newmaxneighbors)
{
	if( newmaxneighbors < currentMaxNeighbors ) return; // If we already have enough memory, there is no need to reallocate...
	// Free the current memory allocation
	delete[] neighborlistspace;

	unsigned neighborlistsize;

	neighborlistsize = newmaxneighbors * (sizeof(int));
	printf("Growing Neighborlist: Atoms = %d  MaxNeighbors/Atom ~ %d  Memory use = %3.2lf Mb \n",
		wspace->atom.size(),
		newmaxneighbors/wspace->atom.size(),
		double(neighborlistsize)/(1024.0*1024.0));

	neighborlistspace = new char[neighborlistsize];
	for(int i = 0; i < wspace->atom.size(); i++) 
	{
		fnbor[i].n = 0;
		fnbor[i].i = (int *) (&neighborlistspace[i * neighborlistsize / wspace->atom.size()]);
		fnbor[i].Type = (int *) (&neighborlistspace[i * neighborlistsize / wspace->atom.size()]);
	}
	currentMaxNeighbors = newmaxneighbors;
	return;
}

void NeighbourList_32Bit_Base::correctOffDiagonalBondOrders()
{
	const BondOrder& bondOrder = wspace->bondorder();
	const std::vector<LongBondOrder>& sparse = bondOrder.getOffDiagData();

	// Loop over all off-diagonal entries
	for( size_t p = 0; p < sparse.size(); p++ )
	{
		int atomIndex = sparse[p].i;
		const std::vector<AtomOrder>& partner = sparse[p].partner;
		// Loop over all its partners
		for( size_t r = 0; r < partner.size(); r++ )
		{
			// Loop over all the neighbors of the atom 
			if( CalcShadow || ( atomIndex > partner[r].j) )
				for(int nj=0; nj < fnbor[atomIndex].n; nj++ )
				{			
					if( NList32Bit_Index(fnbor[atomIndex].i[nj]) == partner[r].j )
					{
						// delete bit field containing the bondorder
						fnbor[atomIndex].i[nj] &= 0xF8FFFFFF;
						// replace those 3 bits with the proper bondorder
						fnbor[atomIndex].i[nj]	|= (partner[r].order&7) << 24;
						break;
					}
				}
			// Loop over all the neighbors of the partner
			if( CalcShadow || ( atomIndex < partner[r].j) )
				for(int nj=0; nj < fnbor[partner[r].j].n; nj++ )
				{			
					if( NList32Bit_Index(fnbor[partner[r].j].i[nj]) == atomIndex )
					{
						// delete bit field containing the bondorder
						fnbor[partner[r].j].i[nj] &= 0xF8FFFFFF;
						// replace those 3 bits with the proper bondorder
						fnbor[partner[r].j].i[nj] |= (partner[r].order&7) << 24;
						break;
					}
				}
		}
	}
}

void NeighbourList::calcNewList()
{
	incFullUpdateCount();   // Mark that we've been called - improtant for other classes to update fully

	// Determine the type of space loaded into workspace to choose
	// the right neighbor list
	if(dynamic_cast<InfiniteSpace *> (&wspace->boundary()) != NULL)
	{
		calcNewList_InfiniteSpace();
	}
	else
		if(dynamic_cast<PeriodicBox *> (&wspace->boundary()) != NULL)
		{
			calcNewList_PeriodicBox();
		}
		else 
		{
			throw(ProcedureException("To use NeighbourList you must load a InfiniteSpace or PeriodicBox space\ninto the workspace you are using, by calling wspace.setBoundary(...) "));
		}

		// finally correct the bondorders that are off the diagonal band matrix
		correctOffDiagonalBondOrders();
}


void NeighbourList::calcNewList_InfiniteSpace()
{
	// Check the right kind of Space is loaded into wspace
	if(dynamic_cast<InfiniteSpace *> (&wspace->boundary()) == NULL){
		throw(ProcedureException("To use NeighbourList_PeriodicBox you must load a PeriodicBox space\ninto the workspace you are using, by calling wspace.setBoundary(...) "));
	}

	// Proxy to atom positions
	const SnapShotAtom *nlist_atom = wspace->cur.atom;

	// Local variables
	int i,j;
	double sqrdistij;
	double sqrIncludeLimit = sqr(Cutoff + Padding);
	dvector dc,ic;

	int bondorderij;
	const BondOrder& bondOrder = wspace->bondorder();

	// restart point for when we run out of neighborlist space
calcNewNeighborList_restart:

	// raw proxies for the bondorder array (faster than std::vector access)
	const int t_MaxIndexDelta = bondOrder.getMaxIndexDelta();
	const char *t_Data  = &bondOrder.getDiagData()[0];

	// clear neighbor arrays
	for(i=0;i<wspace->atom.size();i++) 
		fnbor[i].n=0;

	int memcount = 0;

	// Create neighborlists
	for(i=0;i<wspace->atom.size();i++){				
		//fnbor[i].i = (int *) (&neighborlistspace[memcount*sizeof(int)*2]);
		// set neighbor list start 
		fnbor[i].i = (int *) (&neighborlistspace[memcount*sizeof(int)]);
		fnbor[i].Type = (int *) (&neighborlistspace[memcount*sizeof(int)]);

		int limit;
		if(CalcShadow) limit = wspace->atom.size();
		else           limit = i;

		// now loop over all remaining atoms
		for(j=0;j<limit;j++){
			dc.diff(nlist_atom[j].p,nlist_atom[i].p);
			sqrdistij = dc.innerdot();

			// ignore pairs with distance greater than cutoff
			if(sqrdistij > sqrIncludeLimit) continue;
			// only get bondorder if necessary;
			bondorderij = 7; 
			if( abs(j - i) < t_MaxIndexDelta ){ 
				if( j > i ){
					bondorderij = t_Data[sqrmat(i, j - i - 1, t_MaxIndexDelta)];
				}else{
					bondorderij = t_Data[sqrmat(j, i - j - 1, t_MaxIndexDelta)];
				}
				if( j == i ) continue;
			}

			// when there are more than 12.7 million atoms (unlikely currently)

			//fnbor[i].i[fnbor[i].n*2] = j; // j is a neighbor of i
			//fnbor[i].i[fnbor[i].n*2+1] = (image<<3) + (bondorderij&7);

			// Bit structure:
			// 76543210 76543210 76543210 76543210
			// `---'`-' `-------< index >--------'
			//   |    \------- bondorder            (0-8) 
			//    \----------- image vector (index) (0-32) 

			fnbor[i].i[fnbor[i].n] = j&0x00FFFFFF | ((bondorderij&7)<<24) ;

			fnbor[i].n++;
			memcount++;

		}
		// check if we're starting to run out of space - if so abort, restart and 
		// increase neighborlist
		if(memcount >= (currentMaxNeighbors-wspace->atom.size()-1000) ) 
		{
			reserveMemoryFor( (int)((double)currentMaxNeighbors*1.5));
			goto calcNewNeighborList_restart;
		}
	}
}


void NeighbourList::calcNewList_PeriodicBox()
{
	// Check the right kind of Space is loaded into wspace
	PeriodicBox *periodic_box = dynamic_cast<PeriodicBox *> (&wspace->boundary());
	if(periodic_box == NULL){
		throw(ProcedureException("To use NeighbourList_PeriodicBox you must load a PeriodicBox space\ninto the workspace you are using, by calling wspace.setBoundary(...) "));
	}

	// Proxy to atom positions
	const SnapShotAtom *nlist_atom = wspace->cur.atom;

	// Local variables
	int i,j;
	double sqrdistij;
	double sqrIncludeLimit = sqr(Cutoff + Padding);
	dvector dc,ic;

	int bondorderij;
	const BondOrder& bondOrder = wspace->bondorder();

	int 	image,nimages=wspace->boundary().ncells();
	// restart point for when we run out of neighborlist space
calcNewNeighborList_restart:

	// raw proxies for the bondorder array (faster than std::vector access)
	const int t_MaxIndexDelta = bondOrder.getMaxIndexDelta();
	const char *t_Data  = &bondOrder.getDiagData()[0];

	// clear neighbor arrays
	for(i=0;i<wspace->atom.size();i++) fnbor[i].n=0;

	int memcount = 0;


	// Create neighborlists
	for(i=0;i<wspace->atom.size();i++){				
		//if(i<200) printf("%d  %d\n ",i,memcount);	
		//fnbor[i].i = (int *) (&neighborlistspace[memcount*sizeof(int)*2]);

		// set neighbor list start 
		fnbor[i].i = (int *) (&neighborlistspace[memcount*sizeof(int)]);
		fnbor[i].Type = (int *) (&neighborlistspace[memcount*sizeof(int)]);

		int limit;
		if(CalcShadow) limit = wspace->atom.size();
		else           limit = i;

		// now loop over all remaining atoms
		for(j=0;j<limit;j++){
			dc.diff(nlist_atom[j].p,nlist_atom[i].p);
			image=13; // in Periodic space, image 13 is the self list
			while(dc.x >  periodic_box->halfBoxSize.x){ dc.x -= periodic_box->boxSize.x; image-=9;} 
			while(dc.x < -periodic_box->halfBoxSize.x){ dc.x += periodic_box->boxSize.x; image+=9;} 
			while(dc.y >  periodic_box->halfBoxSize.y){ dc.y -= periodic_box->boxSize.y; image-=3;} 
			while(dc.y < -periodic_box->halfBoxSize.y){ dc.y += periodic_box->boxSize.y; image+=3;} 
			while(dc.z >  periodic_box->halfBoxSize.z){ dc.z -= periodic_box->boxSize.z; image-=1;} 
			while(dc.z < -periodic_box->halfBoxSize.z){ dc.z += periodic_box->boxSize.z; image+=1;}			
			sqrdistij = dc.innerdot();

			// ignore pairs with distance greater than cutoff
			if(sqrdistij > sqrIncludeLimit) continue;

			// only get bondorder if necessary;
			bondorderij = 15; 
			if( abs(j - i) < t_MaxIndexDelta ){ 
				if( j > i ){
					bondorderij = t_Data[sqrmat(i, j - i - 1, t_MaxIndexDelta)];
				}else{
					bondorderij = t_Data[sqrmat(j, i - j - 1, t_MaxIndexDelta)];
				}
				if( j == i ) continue;
				//if(bondorderij < 3) continue;
			}

			//if((i>29)&&(i<31)) printf("%d  %d %d\n ",i,j,memcount);	
			// when there are more than 12.7 million atoms (unlikely currently)

			//fnbor[i].i[fnbor[i].n*2] = j; // j is a neighbor of i
			//fnbor[i].i[fnbor[i].n*2+1] = (image<<3) + (bondorderij&7);

			// Bit structure:
			// 76543210 76543210 76543210 76543210
			// `---'`-' `-------< index >--------'
			//   |    \------- bondorder            (0-8) 
			//    \----------- image vector (index) (0-32) 

			fnbor[i].i[fnbor[i].n] = j&0x00FFFFFF | ((bondorderij&7)<<24) | ((image&31)<<27);

			fnbor[i].n++;
			memcount++;

		}
		// check if we're starting to run out of space - if so abort, restart and 
		// increase neighborlist
		if(memcount >= (currentMaxNeighbors-wspace->atom.size()-1000) ) 
		{
			reserveMemoryFor( (int)((double)currentMaxNeighbors*1.5));
			goto calcNewNeighborList_restart;
		}
	}

	//printf("MEMCOUNT: %d \n", memcount ) ;

}







/// Group based neighbor list
void NeighbourList_GroupBased::calcNewList_PeriodicBox()
{
	// Initial Checks ------------------------------------------------------
	// Check the right kind of Space is loaded into wspace
	PeriodicBox *periodic_box = dynamic_cast<PeriodicBox *> (&wspace->boundary());
	if(periodic_box == NULL){
		throw(ProcedureException("To use NeighbourList_PeriodicBox you must load a PeriodicBox space\ninto the workspace you are using, by calling wspace.setBoundary(...) "));
	}

	// we dont support CalcShadow
	if(CalcShadow){ throw(ProcedureException("This Neighborlist does not support shadow neighbors")); }


	// Actual calculation --------------------------------------------------
	// Proxy to atom positions
	const   SnapShotAtom *nlist_atom = wspace->cur.atom;

	// Local variables
	int     i,j;
	double  sqrdistij;
	double  sqrIncludeLimit = sqr(Cutoff + Padding);
	dvector dc,ic;

	int 	image,nimages=wspace->boundary().ncells();

	// restart point for when we run out of neighborlist space
calcNewNeighborList_restart:

	int   bondorderij;
	const BondOrder& bondOrder = wspace->bondorder();
	const int t_MaxIndexDelta = bondOrder.getMaxIndexDelta();
	const char *t_Data  = &bondOrder.getDiagData()[0];


	int memcount = 0;
	int gi,gj;
	int pairlimit=700;
	int save_group_diff=50;

	// clear neighbor arrays
	for(i=0;i<wspace->atom.size();i++) fnbor[i].n=0;

	for(i=0;i<wspace->atom.size();i++){				
		//if(i<200) printf("%d  %d\n ",i,memcount);	
		fnbor[i].i = (int *) (&neighborlistspace[memcount*sizeof(int)]);
		fnbor[i].Type = (int *) (&neighborlistspace[memcount*sizeof(int)]);
		gi = wspace->atom[i].igroup;
		for(gj=0; gj <=  gi; gj++ ){
			j = wspace->group_index[gj];
			dc.diff(nlist_atom[j].p,nlist_atom[i].p);
			image=13;
			while(dc.x >  periodic_box->halfBoxSize.x){ dc.x -= periodic_box->boxSize.x; image-=9;} 
			while(dc.x < -periodic_box->halfBoxSize.x){ dc.x += periodic_box->boxSize.x; image+=9;} 
			while(dc.y >  periodic_box->halfBoxSize.y){ dc.y -= periodic_box->boxSize.y; image-=3;} 
			while(dc.y < -periodic_box->halfBoxSize.y){ dc.y += periodic_box->boxSize.y; image+=3;} 
			while(dc.z >  periodic_box->halfBoxSize.z){ dc.z -= periodic_box->boxSize.z; image-=1;} 
			while(dc.z < -periodic_box->halfBoxSize.z){ dc.z += periodic_box->boxSize.z; image+=1;}			
			sqrdistij = dc.innerdot();

			if(sqrdistij > sqrIncludeLimit) continue;

			int nnj = wspace->atom.size();
			if(gj<wspace->group_index.size()-1) nnj = wspace->group_index[gj+1];
			for(;j<nnj;j++){
				bondorderij = 7;
				if( abs(j - i) < t_MaxIndexDelta ){ 
					if( j > i ){
						bondorderij = t_Data[sqrmat(i, j - i - 1, t_MaxIndexDelta)];
					}else{
						bondorderij = t_Data[sqrmat(j, i - j - 1, t_MaxIndexDelta)];
					}
					if( j == i ) continue;
					//if(bondorderij < 3) continue;
				}
				//if((i>29)&&(i<31)) printf("%d  %d %d  %d\n ",i,j,memcount,);	
				// Bit structure:
				// 76543210 76543210 76543210 76543210
				// `---'`-' `-------< index >--------'
				//   |    \------- bondorder            (0-8) 
				//    \----------- image vector (index) (0-32) 

				fnbor[i].i[fnbor[i].n] = j&0x00FFFFFF | ((bondorderij&7)<<24) | ((image&31)<<27);
				fnbor[i].n++;
				memcount++;
			}	
		}

		if(memcount >= (currentMaxNeighbors-wspace->atom.size()-1000*30) ) 
		{
			reserveMemoryFor( (int)((double)currentMaxNeighbors*1.5));
			goto calcNewNeighborList_restart;
		}
		continue;

		// now for the remaining atoms in the current group 
		// copy the neighbor list and correct any bondorders..
		int nni = wspace->atom.size();
		if(gi<(wspace->group_index.size()-1)) nni = wspace->group_index[gi+1];
		int isave=1;
		i++;			
		for(;i<nni;i++){
			// make duplicates of the last neighbor list
			fnbor[i].i = (int *) (&neighborlistspace[memcount*sizeof(int)]);
			fnbor[i].Type = (int *) (&neighborlistspace[memcount*sizeof(int)]);
			memcpy( fnbor[i].i, fnbor[i-1].i, fnbor[i-1].n*sizeof(int) );
			fnbor[i].n = fnbor[i-1].n;
			memcount+=fnbor[i].n;
		}
		i = isave;
		for(;i<nni;i++){
			// correct bondorders and other crap 
			for(int nj=0;nj<(min(t_MaxIndexDelta*3,fnbor[i].n));nj++){
				j = fnbor[i].i[nj] & 0x00FFFFFF;
				// turn off inverse neighbors
				if(j>=i){
					fnbor[i].i[nj] &= 0xF8FFFFFF; // set bondor to 0
					continue;
				}
				// correct bondorder for those guys that need it
				if(i - j <= t_MaxIndexDelta ){ 
					bondorderij = t_Data[sqrmat(i, j - i - 1, t_MaxIndexDelta)];
					fnbor[i].i[nj] &= 0xF8FFFFFF; // set bondor to 0
					fnbor[i].i[nj] |= ( (bondorderij&7) << 24 ); // stick in the correct bondorder
				}


			}
		}
		i--;
		// check if we're starting to run out of space - if so abort, restart and 
		// increase neighborlist
		if(memcount >= (currentMaxNeighbors-wspace->atom.size()-1000*30) ) 
		{
			reserveMemoryFor( (int)((double)currentMaxNeighbors*1.5));
			goto calcNewNeighborList_restart;
		}
	}
	//printf("MEMCOUNT: %d \n", memcount ) ;
}




//////// DEPRECATED STUFF /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



/// Vaccum Neighbor List
void DeprecatedNeighbourList::calcNewList()
{
	// check if we're switched on !?
	if(!Enabled) return;
	incFullUpdateCount();   // Mark that we've been called - improtant for other classes to update fully

	// check if we're compatible with the space loaded
	if(dynamic_cast<const InfiniteSpace*>(&wspace->boundary()) == NULL)
	{
		printf("ERROR: DeprecatedNeighbourList is incompatible with the Type of \nERROR: space loaded into workspace (%s)\n", typeid(wspace->boundary()).name());
		throw(ProcedureException("DeprecatedNeighbourList incompatible with the Type of space loaded into workspace"));
	}

	// Proxy to atom positions
	const SnapShotAtom *nlist_atom = wspace->cur.atom;

	// Local variables
	int i,j;
	double sqrdistij;
	double sqrIncludeLimit = sqr(Cutoff + Padding);

	double x,y,z;
	int bondorderij;
	const BondOrder& bondOrder = wspace->bondorder();

	// restart point for when we run out of neighborlist space
calcNewNeighborList_restart:
	// clear neighbor arrays
	for(i=0;i<wspace->atom.size();i++) fnbor[i].n=0;

	// create neighbor lists
	for(i=0;i<wspace->atom.size();i++)
	{
		x = nlist_atom[i].p.x;
		y = nlist_atom[i].p.y;
		z = nlist_atom[i].p.z;

		for(j = i+1;j<wspace->atom.size();j++)
		{
			sqrdistij = sqr((x - nlist_atom[j].p.x)) +
				sqr((y - nlist_atom[j].p.y)) +
				sqr((z - nlist_atom[j].p.z)) ;

			if(sqrdistij < sqrIncludeLimit)
			{
				// only get bondorder if doublely necessary;
				if(!bondOrder.isOffDiag(i,j))
				{
					bondorderij = bondOrder.getDiagBondOrder(i,j);
				}
				else
				{
					bondorderij =  bondOrder.getMaxBondOrder(); // i.e. just normally non-bonded
				}

				// in any case add to full neighbor list
				if(fnbor[i].n >= currentMaxNeighbors) 
				{
					reserveMemoryFor(min((int)wspace->atom.size(), (int)((double)currentMaxNeighbors*1.5)));
					goto calcNewNeighborList_restart;
				}
				fnbor[i].i[fnbor[i].n] = j; // j is a neighbor of i
				if(bondorderij > 3) fnbor[i].Type[fnbor[i].n] = 0;
				else fnbor[i].Type[fnbor[i].n] = 4 - bondorderij;
				fnbor[i].n++;
			}
		}
	}
	// Now analyse the sparse matrix
	const std::vector<LongBondOrder>& sparse = bondOrder.getOffDiagData();
	for( size_t p = 0; p < sparse.size(); p++ )
	{
		int atomIndex = sparse[p].i;
		const std::vector<AtomOrder>& partner = sparse[p].partner;
		for( size_t r = 0; r < partner.size(); r++ )
		{
			for( int q = 0; q < fnbor[atomIndex].n; q++ )
			{			
				if( fnbor[atomIndex].i[q] == partner[r].j )
				{
					fnbor[atomIndex].Type[q] = 4 - partner[r].order;
					break;
				}
			}
		}
		continue;
	}

	if(CalcShadow)
	{
		int fnborsize;
		int nj,itype,jsize;
		for(i=0;i<wspace->atom.size();i++)
		{
			fnborsize = fnbor[i].n;
			for(nj=0;nj<fnborsize;nj++)
			{
				itype = fnbor[i].Type[nj];
				if(itype >= 128) break; // finish if we reached a shadow neighbor
				j = fnbor[i].i[nj];
				if(fnbor[j].n >= currentMaxNeighbors) 
				{
					reserveMemoryFor(min((int)wspace->atom.size(),(int)((double)currentMaxNeighbors*1.5)));
					goto calcNewNeighborList_restart;
				}
				jsize = fnbor[j].n;
				fnbor[j].i[jsize] = i;
				fnbor[j].Type[jsize] = itype + 128;
				fnbor[j].n++;
			}
		}
	}
}


/// Neighbour list that will work for any boundary but not nessessarily particularily
/// fast
void NeighbourList_GeneralBoundary::calcNewList()
{
	// check if we're switched on !?
	if(!Enabled) return;
	incFullUpdateCount();   // Mark that we've been called - improtant for other classes to update fully

	// Proxy to atom positions
	const SnapShotAtom *nlist_atom = wspace->cur.atom;

	// Local variables
	int i,j;
	double sqrdistij;
	double sqrIncludeLimit = sqr(Cutoff + Padding);
	dvector dc;

	int bondorderij;
	const BondOrder& bondOrder = wspace->bondorder();

	// restart point for when we run out of neighborlist space
calcNewNeighborList_restart:

	// clear neighbor arrays
	for(i=0;i<wspace->atom.size();i++) fnbor[i].n=0;


	// Create neighborlists
	for(i=0;i<wspace->atom.size();i++)
	{				
		for(j = i+1;j<wspace->atom.size();j++)
		{
			dc.diff(nlist_atom[i].p,nlist_atom[j].p);
			wspace->boundary().getClosestImage(dc); // only check closest image
			sqrdistij = dc.innerdot();

			if(sqrdistij < sqrIncludeLimit)
			{
				if(!bondOrder.isOffDiag(i,j))
				{ 
					bondorderij = bondOrder.getDiagBondOrder(i,j);
				}
				else
				{
					bondorderij =  bondOrder.getMaxBondOrder(); // i.e. just normally non-bonded
				}
				// in any case add to full neighbor list
				//Save neighbor
				if(fnbor[i].n >= currentMaxNeighbors) 
				{
					reserveMemoryFor(min((int)wspace->atom.size(), (int)((double)currentMaxNeighbors*1.5)));
					goto calcNewNeighborList_restart;
				}

				fnbor[i].i[fnbor[i].n] = j; // j is a neighbor of i
				if(bondorderij > 3) fnbor[i].Type[fnbor[i].n] = 0;
				else fnbor[i].Type[fnbor[i].n] = 4 - bondorderij;
				fnbor[i].n++;
			}
		}
	}


	// Now analyse the sparse matrix
	const std::vector<LongBondOrder>& sparse = bondOrder.getOffDiagData();
	for( size_t p = 0; p < sparse.size(); p++ )
	{
		int atomIndex = sparse[p].i;
		const std::vector<AtomOrder>& partner = sparse[p].partner;
		for( size_t r = 0; r < partner.size(); r++ )
		{
			for( int q = 0; q < fnbor[atomIndex].n; q++ )
			{			
				if( fnbor[atomIndex].i[q] == partner[r].j )
				{
					fnbor[atomIndex].Type[q] = 4 - partner[r].order;
					break;
				}
			}
		}
		continue;
	}


	if(CalcShadow){
		int fnborsize;
		int nj,itype,jsize;
		for(i=0;i<wspace->atom.size();i++){
			fnborsize = fnbor[i].n;
			for(nj=0;nj<fnborsize;nj++){
				itype = fnbor[i].Type[nj];
				if(itype >= 128) break; // finish if we reached a shadow neighbor
				j = fnbor[i].i[nj];
				if(fnbor[j].n >= currentMaxNeighbors) {
					reserveMemoryFor(min((int)wspace->atom.size(),(int)((double)currentMaxNeighbors*1.5)));
					goto calcNewNeighborList_restart;
				}
				jsize = fnbor[j].n;
				fnbor[j].i[jsize] = i;
				fnbor[j].Type[jsize] = itype + 128;
				fnbor[j].n++;
			}
		}
	}


}

/// Neighbour list that will work for any boundary but not nessessarily particularily
/// fast
void NeighbourList_PeriodicBox::calcNewList()
{
	// Check the right kind of Space is loaded into wspace
	PeriodicBox *periodic_box = dynamic_cast<PeriodicBox *> (&wspace->boundary());
	if(periodic_box == NULL){
		throw(ProcedureException("To use NeighbourList_PeriodicBox you must load a PeriodicBox space\ninto the workspace you are using, by calling wspace.setBoundary(...) "));
	}

	// check if we're switched on !?
	if(!Enabled) return;
	incFullUpdateCount();   // Mark that we've been called - improtant for other classes to update fully

	// Proxy to atom positions
	const SnapShotAtom *nlist_atom = wspace->cur.atom;

	// Local variables
	int i,j;
	double sqrdistij;
	double sqrIncludeLimit = sqr(Cutoff + Padding);
	dvector dc;

	int bondorderij;
	const BondOrder& bondOrder = wspace->bondorder();

	// restart point for when we run out of neighborlist space
calcNewNeighborList_restart:

	// clear neighbor arrays
	for(i=0;i<wspace->atom.size();i++) fnbor[i].n=0;

	// Create neighborlists
	for(i=0;i<wspace->atom.size();i++){				
		for(j = i+1;j<wspace->atom.size();j++){
			dc.diff(nlist_atom[i].p,nlist_atom[j].p);
			while(dc.x >  periodic_box->halfBoxSize.x) dc.x -= periodic_box->boxSize.x; 
			while(dc.x < -periodic_box->halfBoxSize.x) dc.x += periodic_box->boxSize.x; 
			while(dc.y >  periodic_box->halfBoxSize.y) dc.y -= periodic_box->boxSize.y; 
			while(dc.y < -periodic_box->halfBoxSize.y) dc.y += periodic_box->boxSize.y; 
			while(dc.z >  periodic_box->halfBoxSize.z) dc.z -= periodic_box->boxSize.z; 
			while(dc.z < -periodic_box->halfBoxSize.z) dc.z += periodic_box->boxSize.z;			
			sqrdistij = dc.innerdot();

			if(sqrdistij < sqrIncludeLimit)
			{
				// only get bondorder if doublely necessary;
				if(!bondOrder.isOffDiag(i,j))
				{
					bondorderij = bondOrder.getDiagBondOrder(i,j);
				}
				else
				{
					bondorderij =  bondOrder.getMaxBondOrder(); // i.e. just normally non-bonded
				}

				// in any case add to full neighbor list
				//Save neighbor
				if(fnbor[i].n >= currentMaxNeighbors) 
				{
					reserveMemoryFor(min((int)wspace->atom.size(), (int)((double)currentMaxNeighbors*1.5)));
					goto calcNewNeighborList_restart;
				}

				fnbor[i].i[fnbor[i].n] = j; // j is a neighbor of i
				if(bondorderij > 3) fnbor[i].Type[fnbor[i].n] = 0;
				else fnbor[i].Type[fnbor[i].n] = 4 - bondorderij;
				fnbor[i].n++;
			}
		}
	}


	// Now analyse the sparse matrix
	const std::vector<LongBondOrder>& sparse = bondOrder.getOffDiagData();
	for( size_t p = 0; p < sparse.size(); p++ )
	{
		int atomIndex = sparse[p].i;
		const std::vector<AtomOrder>& partner = sparse[p].partner;
		for( size_t r = 0; r < partner.size(); r++ )
		{
			for( int q = 0; q < fnbor[atomIndex].n; q++ )
			{			
				if( fnbor[atomIndex].i[q] == partner[r].j )
				{
					fnbor[atomIndex].Type[q] = 4 - partner[r].order;
					break;
				}
			}
		}
		continue;
	}

	if(CalcShadow){
		int fnborsize;
		int nj,itype,jsize;
		for(i=0;i<wspace->atom.size();i++){
			fnborsize = fnbor[i].n;
			for(nj=0;nj<fnborsize;nj++){
				itype = fnbor[i].Type[nj];
				if(itype >= 128) break; // finish if we reached a shadow neighbor
				j = fnbor[i].i[nj];
				if(fnbor[j].n >= currentMaxNeighbors) {
					reserveMemoryFor(min((int)wspace->atom.size(),(int)((double)currentMaxNeighbors*1.5)));
					goto calcNewNeighborList_restart;
				}
				jsize = fnbor[j].n;
				fnbor[j].i[jsize] = i;
				fnbor[j].Type[jsize] = itype + 128;
				fnbor[j].n++;
			}
		}
	}
}


