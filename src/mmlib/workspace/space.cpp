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

#include "system/system.h"

#include "workspace/space.h"

using namespace Maths;


/// makes a box big enough to accomodate the system, adding on a margin
PeriodicBox::PeriodicBox(System &sys, double margin, bool transformToOptimum):
		boxSize(Maths::dvector(1,1,1))  // set to dummy box to start with
{
	// rotate/translate system so as to minimisa box size 
	printf("Creating Auto Periodic Box: Margin: %6.1f , PreOptimise System: %s \n", margin, transformToOptimum ? "true" : "false" );
	if(transformToOptimum){
		
	  sys.zeroCentreOfGeometry();
		sys.alignAlongPrincipalAxes();

	}

	// gets maximum absolute coordinate in system (i.e. largest value of |x|,|y| and |z|
	boxSize = sys.getEncompassingVector();  
	boxSize.mul(2);
	boxSize.add( dvector(margin,margin,margin) );
	init();
}

void PeriodicBox::setupCellDimensions(double Cutoff)
{
	double maxradius;   
	
	m_cutoff = Cutoff;
	m_sqrcutoff = sqr(Cutoff);

	// maximal radius of sphere that fis completely into box
	dvector limit(halfBoxSize);
	maxradius = min(limit.x,limit.y,limit.z);

	// decide if we need replicas of boxes ?
	need_replicas = Cutoff < maxradius ? false : true;

	nx = (unsigned)ceil(Cutoff/limit.x);
	ny = (unsigned)ceil(Cutoff/limit.y);
	nz = (unsigned)ceil(Cutoff/limit.z);
	
	cellSize.setTo(
		double(nx)*boxSize.x,
		double(ny)*boxSize.y,
		double(nz)*boxSize.z
	);
	halfCellSize.setTo(cellSize);
	halfCellSize.div(2.0);

 	unsigned nreplicas = nx * ny * nz ;  //-1 ?

	celloffset.resize(nreplicas);
	
	int i,j,k,count_image=0;

	for(i=0;i<nx;i++){
		for(j=0;j<ny;j++){
			for(k=0;k<nz;k++){
				celloffset[count_image].setTo(
					double(i)*boxSize.x,
					double(j)*boxSize.y,
					double(k)*boxSize.z
				);
				count_image++;
			}
		}
	}

	basisvector.resize( (2*nx+1)*(2*ny+1)*(2*ny+1) );
	count_image=0;
	for(int i=-1;i<=int(nx);i++){
		for(int j=-1;j<=int(ny);j++){
			for(int k=-1;k<=int(nz);k++){
				basisvector[count_image].setTo(
					double(i)*boxSize.x,
					double(j)*boxSize.y,
					double(k)*boxSize.z
				);
				//printf("%f  %f  %f \n", double(i)*boxSize.x,
				//         double(j)*boxSize.y,
				//					          double(k)*boxSize.z
				//V										);
				count_image++;
			}
		}
	}
}


void PeriodicBox::info() const{
	printf("Space: Periodic Box  \n  Dimensions: <%5.1f,%5.1f,%5.1f> \n  MaxImages:  <%d,%d,%d>  \n  Cell:      <%5.1f,%5.1f,%5.1f>\n",
		boxSize.x,
		boxSize.y,
		boxSize.z,
		nx,ny,nz,
		cellSize.x,
		cellSize.y,
		cellSize.z);				
	printf("  Volume:   %7.0f A^3 \n", volume() );	
}

// provides a set of points even distributed through the space
void PeriodicBox::getEvenPointDistribution(
	std::vector<dvector> &pointarray, 
	unsigned N
	) const 
{
	pointarray.clear();
	unsigned i,j,k;
	// factorise N into three numbers

	double dNx = pow( sqr(boxSize.x) / ( boxSize.y * boxSize.z ) * double(N) , 1.0/3.0);
	double dNy = pow( sqr(boxSize.y) / ( boxSize.x * boxSize.z ) * double(N) , 1.0/3.0);
	double dNz = pow( sqr(boxSize.z) / ( boxSize.x * boxSize.y ) * double(N) , 1.0/3.0);

	//printf(" %e %e %e \n",dNx,dNy,dNz);

	// now find closest whole three numbers to the required N

	unsigned iNx = (unsigned)ceil(dNx);
	unsigned iNy = (unsigned)ceil(dNy);
	unsigned iNz = (unsigned)ceil(dNz);

	std::vector<double> Nx_choices; 
	std::vector<double> Ny_choices; 
	std::vector<double> Nz_choices; 

	unsigned bestNx=iNx;
	unsigned bestNy=iNy;
	unsigned bestNz=iNz;
	unsigned best = bestNx*bestNy*bestNz - N;   // the excess is to be minimised.
	for(i=iNx-2;i<=iNx+4;i++){
		for(j=iNy-2;j<=iNy+4;j++){
			for(k=iNz-2;k<=iNz+4;k++){
				unsigned tbest = i*j*k - N;
				if((tbest<best)){
					best = tbest;
					bestNx=i;
					bestNy=j;
					bestNz=k;
				}
			}
		}
	}

	//printf(" %d %d %d \n",bestNx,bestNy,bestNz);

	dvector Step(boxSize);
	Step.x /= double(bestNx);
	Step.y /= double(bestNy);
	Step.z /= double(bestNz);

	for(i=0;i<bestNx;i++){
		for(j=0;j<bestNy;j++){
			for(k=0;k<bestNz;k++){
				dvector newpoint(Step.x*(double(i)+0.5),Step.y*(double(j)+0.5),Step.z*(double(k)+0.5));
				newpoint.sub(halfBoxSize);
				pointarray.push_back(newpoint);
				if(pointarray.size() >= N) return;
			}
		}
	}
}


ClosedSpace *castTo_ClosedSpace( Space *sp ){
	return dynamic_cast<ClosedSpace *> (sp);
}


