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

#include "umbrella.h"

namespace Monitors{

	void printUmbrellaProfile(UmbrellaMonitor &umbrella, double binwidth){

		Maths::Histogram1D hist = umbrella.getHistogram(binwidth);
		std::deque<double> probability;
		hist.getProb(probability);

		for(size_t x = 0; x < hist.data_size(); x++) {
			double Q = hist.x1() + (double) x * (double) hist.xgridsize();

			printf("%shst%10.5lf\t%8d\t%10.8lf\t%10.8lf\t%10.8lf\n",
				umbrella.name.c_str(),
				Q,
				int(hist.data(x)),
				probability[x],
				-0.592*log(probability[x]),
				-0.592*log(probability[x]) - umbrella.ff->calcEnergyAtQ(Q) );
		}
	}
}

