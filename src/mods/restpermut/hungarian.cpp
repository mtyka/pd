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

#include <valarray>

// used to 'simulate' 2D arrays using linear arrays
#define sqrmati(x,y,a) ((y)*(a) + (x))
#include "exception.h"
#include "maths/maths.h"

// Hungarian Algorithm for Linear assignement - the original 
// algorithm due to Kuhn and Munkres. Its pretty slow ( O(n^4) )
// 
// Kuhn, H. W.: The Hungarian method for the assignment problem. 
// Naval Research Logistics Quarterly 2, 83-97 (1955).

int LinearAssignmentHungarian(
	std::valarray<int> &costm, 
	std::vector<int> &row2col, 
	std::vector<int> &col2row
){
	using namespace Maths;
	if( row2col.size() != col2row.size() ){
		throw(CodeException("LinearAssignmentHungarian: row2col and col2row must be of equal size"));
	}
	if( sqr(row2col.size()) != costm.size() ){
		throw(CodeException("LinearAssignmentHungarian: row2col size squared must be equal to size of costmatrix"));
	}

  
  int           i,j,irow,icol,irowstart;
  unsigned char *smat = new unsigned char[costm.size()];
  unsigned char *colmark = new unsigned char[row2col.size()];
  unsigned char *rowmark = new unsigned char[row2col.size()];
  int           ncovcols;
  int           *Srow = new int[row2col.size()];
  int           *Scol = new int[row2col.size()];
  int           *Prow = new int[row2col.size()];
  int           *Pcol = new int[row2col.size()];
  int           iS,nS;
  int           iP,nP;
  int           low;
  int           ilow;
  int           cnt;
  int           pos;

  int           stat0=0;
  int stat1=0;
  // Step 1 - find smallest value in each row & substract from rest

  const unsigned char zero        = 0x01;
  const unsigned char unmarkedzero = 0x02;
  const unsigned char starred      = 0x04;
  const unsigned char primed       = 0x08;
  const unsigned char rowline    = 0x10;
  const unsigned char colline    = 0x20;
  const unsigned char unused1    = 0x40;
  const unsigned char unused2    = 0x80;

	size_t nlist = row2col.size();
  for(i=0;i<(nlist*nlist);i++)smat[i]=0;
  for(i=0;i<nlist;i++)colmark[i]=0;
  for(i=0;i<nlist;i++)rowmark[i]=0;

  for(irow=0;irow<nlist;irow++){
    //find lowest value in each row and substract it from all the values of that row
    irowstart=sqrmati(irow,0,nlist);
    low = costm[sqrmati(irow,0,nlist)];
    ilow = 0;
    for(icol=1;icol<nlist;icol++){
      if( costm[sqrmati(irow,icol,nlist)] < low){
        low = costm[sqrmati(irow,icol,nlist)];
        ilow = icol;
      } 
    }
    for(icol=0;icol<nlist;icol++) costm[sqrmati(irow,icol,nlist)] -= low;
  
    smat[sqrmati(irow,ilow,nlist)] |= zero; // mark the zero  
    // Step 2
    // also mark the zero (each row should contain one) unless there is already a zero in that row
    // i.e, first check if there's a zero already;
    for(i=0;i<irow;i++){ // for every row thus far handled
     if( (smat[sqrmati(i,ilow,nlist)] & starred) != 0 ) break;
    }
    if(i>=irow) smat[sqrmati(irow,ilow,nlist)] |= starred; // if no star found above, star this zero
  }


  // Step 3 - mark all columns that contain starred zeros with a line
  // if all columns contain a star, this is the optimum assignment and algoritm will quit

  do{
    
    ncovcols = 0;
 
    for(icol=0;icol<nlist;icol++)
		{
      int foundstar=0;
      for(irow=0;irow<nlist;irow++)
			{
				if( (smat[sqrmati(irow,icol,nlist)] & starred) != 0 ) 
				{ 
					foundstar = 1; 
					break; 
				} 
      }
      if(foundstar==0) continue;
      colmark[icol] = 1;     
      ncovcols++;
    } 

    if(ncovcols == nlist) break; // finish if al lcolumns are covered
   
    // Step 4

    step4:
    ilow=-1;
    low = 100000000;
    icol=0;
    irow=0;
    nP=0;
    nS=0;

    pos=0;
    for(icol=0;icol<nlist;icol++)
		{
     if(colmark[icol] != 0)
		 { 
			 pos+= nlist; // skip marked columns
			 continue; 
		 }   
     for(irow=0;irow<nlist;irow++)
		 {
       if(rowmark[irow] != 0)
			 { 
				 pos++; // skip marked rows
				 continue; 
			 }  
       if((costm[pos]<low)){   // record lowest incase we dont find a zero
         low = costm[pos];
         ilow = pos;
       }
       
       if( (smat[pos] & zero) != 0 ){    // if find an uncovered zero
         smat[pos] |= primed;            // prime it 
         for(i=0;i<nlist;i++){                               // search the row for a starred 0
           if( (smat[sqrmati(irow,i,nlist)] & starred ) != 0){// if find starred 0
               rowmark[irow]=1;                              // cover it's row
               colmark[i]=0;                                 // uncover it's column
               goto step4;                                   // repeat Step 4 (goto ok here!)
           }
         }                                                   // if no starred zero in row
         Prow[nP] = irow;                                    // remember row&column as P0
         Pcol[nP] = icol;
         nP++;
         break;
       }
       pos++;
     }
     if(nP>0) break;
    } 

    if(nP==0){
      stat1++;
      // cannot find an uncovered zero so do part 6
      smat[ilow] |= zero;                                      // mark new zero (to be created by subsequent subtraction) 
      
      unsigned tpos=0;
			unsigned nlistd4 = nlist/4;
			unsigned nlistm4 = nlist%4;
			unsigned nlistcnt;
			for(icol=0;icol<nlist;icol++)
			{
				if(colmark[icol] != 0){
					// hand-unrolled loop
					//for(irow=0;irow<nlist;irow++) 
					irow=0;
					for(nlistcnt=0;nlistcnt<nlistd4;nlistcnt++)
					{
						if(rowmark[irow] != 0)						{ 
							costm[tpos] += low;             // add lowest uncovered value to doubly covered squares
							smat[tpos] &= ~zero;           // it consequently cannot be 0 anymore 
						}
						tpos++; irow++;
						if(rowmark[irow] != 0)						{ 
							costm[tpos] += low;             // add lowest uncovered value to doubly covered squares
							smat[tpos] &= ~zero;           // it consequently cannot be 0 anymore 
						}
						tpos++; irow++;
						if(rowmark[irow] != 0)						{ 
							costm[tpos] += low;             // add lowest uncovered value to doubly covered squares
							smat[tpos] &= ~zero;           // it consequently cannot be 0 anymore 
						}
						tpos++; irow++;
						if(rowmark[irow] != 0)						{ 
							costm[tpos] += low;             // add lowest uncovered value to doubly covered squares
							smat[tpos] &= ~zero;           // it consequently cannot be 0 anymore 
						}
						tpos++; irow++;
					}
					for(nlistcnt=0;nlistcnt<nlistm4;nlistcnt++){
						if(rowmark[irow] != 0)
						{ 
							costm[tpos] += low;             // add lowest uncovered value to doubly covered squares
							smat[tpos] &= ~zero;           // it consequently cannot be 0 anymore 
						}
						tpos++; irow++;
					}
				}else{
					irow=0;

					for(nlistcnt=0;nlistcnt<nlistd4;nlistcnt++)
					{
						if((rowmark[irow] == 0)){
							costm[tpos] -= low;             // subtract lowest uncovered value from uncovered squares
						}
						tpos++;irow++;
						if((rowmark[irow] == 0)){
							costm[tpos] -= low;             // subtract lowest uncovered value from uncovered squares
						}
						tpos++;irow++;
						if((rowmark[irow] == 0)){
							costm[tpos] -= low;             // subtract lowest uncovered value from uncovered squares
						}
						tpos++;irow++;
						if((rowmark[irow] == 0)){
							costm[tpos] -= low;             // subtract lowest uncovered value from uncovered squares
						}
						tpos++;irow++;
					}
					for(nlistcnt=0;nlistcnt<nlistm4;nlistcnt++){
						if((rowmark[irow] == 0)){
							costm[tpos] -= low;             // subtract lowest uncovered value from uncovered squares
						}
						tpos++;irow++;
					}
					
				}


      }
      goto step4;
    }

    // Step 5
    // print matrix:
    stat0++;
    iS=0;
    iP=0;
    do{
      // find a starred S[iS] in P[iP]'s column
      iS++;
      for(irow=0;irow<nlist;irow++){
        if( (smat[sqrmati(irow,Pcol[iP],nlist)] & starred) != 0 ){    // if find an star
          Srow[iS]=irow;
          Scol[iS]=Pcol[iP];
          break;
        } 
      }
      if(irow>=nlist){                                                // couldn't find a star in P[iP]'s column
         
         break;						   // finish and process series P0,S1,P1,S2,...PN 
      }
      iP++;
      // find a new primed P[iP] in S[iS]'s row
      for(icol=0;icol<nlist;icol++){
        if( (smat[sqrmati(Srow[iS],icol,nlist)] & primed) != 0 ){    // if find an prime
          Prow[iP]=Srow[iS];
          Pcol[iP]=icol;
          break;
        } 
      }
      if(icol>=nlist){                                                // couldn't find a star in P[iP]'s column
          printf("Assert error: This shouldn't occur \n");
          exit(0);
      }
  

    }while(1);

    nS = iS-1;
    nP = iP;
    //  process series P0,S1,P1,S2,...PN 
    for(iS=1;iS<=nS;iS++) smat[sqrmati(Srow[iS],Scol[iS],nlist)] &= ~starred; // unstar each starred zero in series
    //  process series P0,S1,P1,S2,...PN            
    for(iP=0;iP<=nP;iP++) smat[sqrmati(Prow[iP],Pcol[iP],nlist)] |=  starred; // star each primed zero in series
    
    // unprime everything
    for(i=0;i<(nlist*nlist);i++) smat[i] &= ~primed;  
    // uncover all the rows/cols
    for(i=0;i<nlist;i++){
      colmark[i]=0;
      rowmark[i]=0;
    }
  } while(1);

  // now set all the assignments
  cnt=0;
  for(i=0;i<nlist;i++){
   for(j=0;j<nlist;j++){
    if((smat[sqrmati(j,i,nlist)] & starred)!=0){
      row2col[j] = i;
      col2row[i] = j;
      cnt+=nlist-j;
      break;
    }
    cnt++;
   }
  }
  
  //printf("%d  %d \n",stat0,stat1);
  
	// clean up 
	delete [] smat; 
  delete [] colmark;
  delete [] rowmark;
  delete [] Srow; 
  delete [] Scol; 
  delete [] Prow; 
  delete [] Pcol; 

	// work out final cost
	int cost = 0;
	for(i=0; i<nlist; i++)
	{
		j = col2row[i];
		cost += costm[sqrmat(i,j,nlist)];
	}

	return cost;
}


//  Adapted from PASCAL version in
//  Jonker, R. and Volgenant, A., 1987. A shortest augmenting
//  path algorithm for dense and sparse linear assignment problems. 
//  Computing 38, pp. 325–340
int LinearAssignmentJVC(
	std::valarray<int> &costm,
	std::vector<int> &row2col,
	std::vector<int> &col2row
)
{
	// asserts
	if( row2col.size() != col2row.size() ){
		throw(CodeException("LinearAssignmentHungarian: row2col and col2row must be of equal size"));
	}
	if( Maths::sqr(row2col.size()) != costm.size() ){
		throw(CodeException("LinearAssignmentHungarian: row2col size squared must be equal to size of costmatrix"));
	}

	size_t nlist = col2row.size();
	bool found_unassigned;
	int i, i1, f0=0, lastfree, f, i0, k, f1;
	int j, j1, j2, end, last, low, up;
	int min, h, u1, u2, tmp1;
	int *pred    = new int[nlist];	//predecessor-array for shortest path tree
	int *free    = new int[nlist]; //unassigned rows (number f0, index f)
	int *col     = new int[nlist]; //col array of columns, scanned
	int *match	 = new int[nlist];
	int *d       = new int[nlist];      // shortest path lengths
	int *u       = new int[nlist];   // row cost
	int *v       = new int[nlist];   // column cost

	for(i=0; i<nlist; i++) match[i] = 0;

	//column reduction
	for(j=nlist-1; j>=0; j--)
	{
		min = costm[sqrmat(0,j,nlist)];
		i1 = 0;
		for(i=1; i<nlist; i++)
			if(costm[sqrmat(i,j,nlist)]<min)
			{
				min = costm[sqrmat(i,j,nlist)];
				i1 = i;
			}
			v[j] = min;

			if(++match[i1]==1)
			{
				col2row[i1] = j;
				row2col[j] = i1;
			}else{
				row2col[j] = -1;
			}
	}

	// reduction transfer
	for(i=0; i<nlist; i++){
		if(match[i]==0){			
			free[f0++] = i;
		}else{
			if(match[i]==1)
			{
				j1 = col2row[i];
				min = INT_MAX;
				for(j=0; j<nlist; j++)
					if(j!=j1){
						if(costm[sqrmat(i,j,nlist)]-v[j] < min){
							min = costm[sqrmat(i,j,nlist)] - v[j];
						}
					}
				v[j1] = v[j1] - min;
			}
		}
	}

	//augmenting row reduction
	int loopcnt = 0;
	do
	{
		loopcnt++;
		k=0;
		lastfree = f0;
		f0 = 0;
		while(k<lastfree)
		{
			i = free[k];
			k++;

			u1 = costm[sqrmat(i,0,nlist)] - v[0];
			j1 = 0;
			u2 = INT_MAX;
			for(j=1; j<nlist; j++)
			{
				h = costm[sqrmat(i,j,nlist)] - v[j];
				if(h<u2){
					if(h>=u1)
					{
						u2 = h;
						j2 = j;
					}
					else
					{
						u2 = u1;
						u1 = h;
						j2 = j1;
						j1 = j;
					}
				}
			}

			i0 = row2col[j1];
			if(u1 < u2){
				v[j1] = v[j1] - (u2 - u1);
			}else{
				if(i0>=0)
				{
					j1 = j2;
					i0 = row2col[j2];
				}
			}

			col2row[i] = j1;
			row2col[j1] = i;

			if(i0 >= 0){
				if(u1 < u2){
					free[--k] = i0;
				}else{
					free[f0++] = i0;
				}
			}
		}
	}
	while(loopcnt < 2);  // routine applied twice

	//augmentation
	for(f=0; f<f0; f++)
	{
		f1 = free[f];
		low = 0;
		up = 0;

		for(j=0; j<nlist; j++)
		{
			d[j] = costm[sqrmat(f1,j,nlist)] - v[j];
			pred[j] = f1;
			col[j] = j;
		}

		found_unassigned = false;
		do
		{
			if(up==low)
			{
				last = low - 1;

				min = d[col[up++]];
				for(k=up; k<nlist; k++)
				{
					j = col[k];
					h = d[j];
					if(h<=min)
					{
						if(h<min)
						{
							up = low;
							min = h;
						}
						col[k] = col[up];
						col[up++] = j;
					}
				}
				for(k=low; k<up; k++){
					if(row2col[col[k]] < 0)
					{
						end = col[k];
						found_unassigned = true;
						break;
					}
				}
			}

			if(!found_unassigned)
			{
				j1 = col[low];
				low++;
				i = row2col[j1];
				h = costm[sqrmat(i,j1,nlist)] - v[j1] - min;

				for(k=up; k<nlist; k++)
				{
					j = col[k];
					tmp1 = costm[sqrmat(i,j,nlist)] - v[j] - h;
					if(tmp1<d[j])
					{
						pred[j] = i;
						if(tmp1==min)
							if(row2col[j]<0)
							{
								end = j;
								found_unassigned = true;
								break;
							}
							else
							{
								col[k] = col[up];
								col[up++] = j;
							}
						d[j] = tmp1;
					}
				}
			}
		}
		while(!found_unassigned);

		for(k=0; k<=last; k++)
		{
			j1 = col[k];
			v[j1] = v[j1] + d[j1] - min;
		}

		do
		{
			i = pred[end];
			row2col[end] = i;
			j1 = end;
			end = col2row[i];
			col2row[i] = j1;
		}
		while(i!=f1);
	}

	// work out final cost
	int cost = 0;
	for(i=0; i<nlist; i++)
	{
		j = col2row[i];
		u[i] = costm[sqrmat(i,j,nlist)] - v[j];
		cost += costm[sqrmat(i,j,nlist)];
	}

	delete [] pred;
	delete [] free;
	delete [] col;
	delete [] match;
	delete [] d;
	delete [] u;
	delete [] v;

	return cost;
}

