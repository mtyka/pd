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

#include "workspace/workspace.h"

#include "monitors/monitorbase.h"

namespace Monitors{

	MonitorBase::MonitorBase() 
	{
		name = "MonitorBase";
		formatstring = " %10.7e";
		headerformatstring = " %14s";
		OutputLevel = Verbosity::Silent;
		PrintAverage = false;
		granularity = 1;
		curdata=0;
		gran_buffer=0;
		gran_count=0;
		lastdata=0;
		lastdata_valid=false;

		RawIgnoreFirstN = 0;
		RawSkip = 1;
		raw_count = 0;
	}

	int MonitorBase::measure()  
	{
		// Should we record any data at all - the raw_xx feature 
		// allows ignoring the first n measurements and only recording every so and so-th
		// data.
		raw_count++;
		if(((int)raw_count-1)<RawIgnoreFirstN) return 0;
		if(((int)(raw_count-1)%RawSkip)!=0) return 0;

		// call overloaded implementation of the derived class to make the measurement
		setcurdata();                   // make derived class work out it's value
		return 0;
	}

	int MonitorBase::addData(double newdata)  
	{
		curdata = newdata;
		lastdata = curdata;
		lastdata_valid = true;
		if(granularity==1){							// gran. of 1 means each piece of data is stored individually
			dataarray.push_back(curdata);
		}else{													// otherwise the data is 'bundled' in groups of 'granularity'
			gran_buffer += curdata;				// and stored as averages over those bundles. This does not affect
			gran_count++;									// overall averages but does decrease any calculated standard deviation
			if(gran_count >= granularity){ // of the data. SDs over block averages (bootstrapping) are also unaffected.
				dataarray.push_back(gran_buffer/granularity);
				gran_buffer=0;
				gran_count=0;
			}
		}
		return 0;
	}

	// repeat last measurement (if there was one)
	int MonitorBase::repeat()  
	{
		if(!lastdata_valid) return 0;

		// Should we record any data at all - the raw_xx feature 
		// allows ignoring the first n measurements and only recording every so and so-th
		// data.
		raw_count++;
		if(((int)raw_count-1)<RawIgnoreFirstN) return 0;
		if(((int)(raw_count-1)%RawSkip)!=0) return 0;
		//printf("DEBUG LASTDATA: %f \n",lastdata);
		if(granularity==1){							// gran. of 1 means each piece of data is stored individually
			dataarray.push_back(lastdata);
		}else{													// otherwise the data is 'bundled' in groups of 'granularity'
			gran_buffer += lastdata;				// and stored as averages over those bundles. This does not affect
			gran_count++;									// overall averages but does decrease any calculated standard deviation
			if(gran_count >= granularity){ // of the data. SDs over block averages (bootstrapping) are also unaffected.
				dataarray.push_back(gran_buffer/granularity);
				gran_buffer=0;
				gran_count=0;
			}
		}
		return 0;
	}

	void MonitorBase::reset()
	{
		curdata=0;
		gran_buffer=0;
		gran_count=0;
		raw_count = 0;
		dataarray.clear();
	}

	void MonitorBase::printHeader() const
	{
		if(OutputLevel){
			printf((char*)headerformatstring.c_str(),(const char*)name.c_str());
			if(PrintAverage){
				printf((char*)headerformatstring.c_str(),"<-- Av.");
			}
		}
	}

	void MonitorBase::printCurData() const
	{
		if(OutputLevel){
			printf(formatstring.c_str(),getCurData());
			if(PrintAverage){
				printf(formatstring.c_str(),getAverageSubdata() );
			}
		}
	}

	void MonitorBase::printData(size_t index) const
	{
		double pdata;
		if(index>=dataarray.size()) pdata = Maths::dNAN();
		else pdata = dataarray[index];
		printf(formatstring.c_str(),pdata);
	}

	void MonitorBase::printAllData() const
	{
		for(int i=0;i<dataarray.size();i++){
			printf("%s %.7e\n",name.c_str(),dataarray[i]);
		}
	}

	void MonitorBase::printHistogram(double binwidth)
	{
		printHistogram(binwidth, 0, dataarray.size());
	}

	void MonitorBase::printHistogram(double binwidth, unsigned start, unsigned end)const 
	{
		Maths::Histogram1D hist;

		hist = getHistogram(binwidth, start, end);
		hist.print_count();
	}

	double MonitorBase::getData(size_t index)const 
	{
		if(index > dataarray.size()) return Maths::dNAN();
		return dataarray[index];
	}

	double MonitorBase::getAverageSubdata()const 
	{
		return getAverageSubdata(0,dataarray.size());
	}
	double MonitorBase::getAverageSubdata(size_t start)const 
	{
		return getAverageSubdata(start,dataarray.size());
	}

	double MonitorBase::getAverageSubdata(size_t start, size_t end)const 
	{
		double sum = 0;
		for(unsigned i = start;i<end;i++){
			sum += dataarray[i];
		}
		return sum/double(end-start);
	}

	double MonitorBase::getVarianceSubdata(size_t start, size_t end)const 
	{
		double sum = 0;
		double sumsqr = 0;
		for(size_t i = start;i<end;i++){
			sum += dataarray[i];
			sumsqr += Maths::sqr(dataarray[i]);
		}
		return (sumsqr/double(end-start)) - Maths::sqr((sum/double(end-start))) ;
	}


	double MonitorBase::getAverage(double start, double end)const 
	{
		if((start<0.0)||(start>1.0)||
			(end<0.0)||(end>1.0)){
				THROW(ArgumentException,"ERROR: start and end must be between 0.0 and 1.0");
		}

		return getAverageSubdata(unsigned(start*double(nData())),
			unsigned(end*double(nData())));
	}


	double MonitorBase::getVariance(double start, double end)const 
	{
		if((start<0.0)||(start>1.0)||
			(end<0.0)||(end>1.0)){
				THROW(ArgumentException,"ERROR: start and end must be between 0.0 and 1.0");
		}

		return getVarianceSubdata(unsigned(start*double(nData())),
			unsigned(end*double(nData())));
	}

	/// Performs a bootstrap block analysis of the averages etc..
	void MonitorBase::printBlockAverage(double ignore, unsigned blocks, int blocks_end)
	{
		unsigned nblocks;

		argcheck_range("ignore",ignore,0.0,1.0);
		argcheck_gt("blocks",blocks,(unsigned)1);

		if(((int)blocks_end) < ((int)blocks)) blocks_end = blocks;
		if(blocks_end == blocks){
			printf("bs %s av % 8.4lf | ",name.c_str(), getAverage(ignore,1.0) );
		}else{
			printf("bs %s av % 8.4lf\n",name.c_str(), getAverage(ignore,1.0) );
		}

		for(nblocks = blocks; nblocks <= blocks_end; nblocks++){ // for various number of blocks
			if(blocks_end != blocks){
				printf("bs %s %2d ",name.c_str(), nblocks);
			}

			double sumx=0.0;
			double sumx2=0.0;
			unsigned b;
			for(b=0;b<nblocks; b++){ // calculate average over each block
				double bav = getAverage(
					Maths::min(1.0,ignore + double(b)*(1.0 - ignore)/double(nblocks) ) ,
					Maths::min(1.0,ignore + double(b+1)*(1.0 - ignore)/double(nblocks) ) );
				printf("% 8.4lf ", bav );
				sumx += bav;
				sumx2 += Maths::sqr(bav);
			}
			for(;b<blocks_end;b++){
				printf(" ");
			}
			printf("stddev: % 8.4lf", sqrt( (sumx2 - Maths::sqr(sumx)/double(nblocks))/(nblocks-1) ) );
			printf("\n");
		}
	}

	/// Performs a bootstrap block analysis of the averages etc..
	double MonitorBase::getBlockAvStdDev(double ignore, unsigned nblocks)
	{
		argcheck_range("ignore",ignore,0.0,1.0);
		argcheck_gt("blocks",nblocks,(unsigned)1);

		double sumx=0.0;
		double sumx2=0.0;
		unsigned b;
		for(b=0;b<nblocks; b++){ // calculate average over each block
			double bav = getAverage(
				Maths::min(1.0,ignore + double(b)*(1.0 - ignore)/double(nblocks) ) ,
				Maths::min(1.0,ignore + double(b+1)*(1.0 - ignore)/double(nblocks) ) );
			sumx += bav;
			sumx2 += Maths::sqr(bav);
		}
		return sqrt( (sumx2 - Maths::sqr(sumx)/double(nblocks))/(nblocks-1) );
	}

	void MonitorBase::printRunningAverageSubdata(unsigned start, unsigned end, unsigned Step)const 
	{
		for(unsigned finish = (start+Step); finish < end; finish += Step )
		{
			printf("%s %.7e\n",name.c_str(), getAverageSubdata(start,finish) );
		}
	}


	void MonitorBase::printRunningAverage(double start, double end, unsigned Step)const 
	{
		if((start<0.0)||(start>1.0)||
			(end<0.0)||(end>1.0)){
				THROW(ArgumentException,"ERROR: start and end must be between 0.0 and 1.0");
		}

		return printRunningAverageSubdata(unsigned(start*double(dataarray.size())),
			unsigned(end*double(dataarray.size())),
			Step);
	}


	Maths::Histogram1D MonitorBase::getHistogram(
		double binwidth, 
		unsigned start
		)const
	{
		return getHistogram(binwidth, start, dataarray.size());
	}

	Maths::Histogram1D MonitorBase::getHistogram(
		double binwidth,
		unsigned start,
		unsigned end
		)const 
	{
		Maths::Histogram1D hist(binwidth);

		for(unsigned i = start;i<end;i++){
			hist.addPoint(dataarray[i]);
		}
		return hist;
	}

	double MonitorBase::getAverageSubdata(
		const MonitorBase &weight,
		unsigned start,
		unsigned end
		)const 
	{
		double sum = 0;
		double weightsum = 0;
		for(unsigned i = start;i<end;i++){
			sum += dataarray[i]*weight.getData(i);
			weightsum += weight.getData(i);
		}
		return sum/weightsum;
	}

	double MonitorBase::getAverageSubdata(const MonitorBase &weight)const 
	{
		return getAverageSubdata(weight,0,dataarray.size());
	}

	double MonitorBase::getAverageSubdata(const MonitorBase &weight,unsigned start)const
	{
		return getAverageSubdata(weight,start,dataarray.size());
	}

	double MonitorBase::getAverage(const MonitorBase &weight, double start, double end)const 
	{
		if((start<0.0)||(start>1.0)||
			(end<0.0)||(end>1.0)){
				THROW(ArgumentException,"ERROR: start and end must be between 0.0 and 1.0");
		}

		return getAverageSubdata(weight,
			unsigned(start*double(nData())),
			unsigned(end*double(nData())));
	}

	/// Performs a bootstrap block analysis of the averages etc..
	void MonitorBase::printBlockAverage(
		const MonitorBase &weight, 
		double ignore, 
		unsigned blocks, 
		int blocks_end
		)
	{
		unsigned nblocks;

		argcheck_range("ignore",ignore,0.0,1.0);
		argcheck_gt("blocks",blocks,(unsigned)1);

		if(((int)blocks_end) < ((int)blocks)) blocks_end = blocks;
		if(blocks_end == blocks){
			printf("bs %s av % 8.4lf | ",name.c_str(), getAverage(weight, ignore,1.0) );
		}else{
			printf("bs %s av % 8.4lf\n",name.c_str(), getAverage(weight, ignore,1.0) );
		}

		for(nblocks = blocks; nblocks <= blocks_end; nblocks++){ // for various number of blocks
			if(blocks_end != blocks){
				printf("bs %s %2d ",name.c_str(), nblocks);
			}

			double sumx=0.0;
			double sumx2=0.0;
			unsigned b;
			for(b=0;b<nblocks; b++){ // calculate average over each block
				double bav = getAverage(weight,
					Maths::min(1.0,ignore + double(b)*(1.0 - ignore)/double(nblocks) ) ,
					Maths::min(1.0,ignore + double(b+1)*(1.0 - ignore)/double(nblocks) ) );
				printf("% 8.4lf ", bav );
				sumx += bav;
				sumx2 += Maths::sqr(bav);
			}
			for(;b<blocks_end;b++){
				printf(" ");
			}
			printf("stddev: % 8.4lf", sqrt( (sumx2 - Maths::sqr(sumx)/double(nblocks))/(nblocks-1) ) );
			printf("\n");
		}

	}


	/// Performs a bootstrap block analysis of the averages etc..
	double MonitorBase::getBlockAvStdDev(
		const MonitorBase &weight, 
		double ignore, 
		unsigned nblocks
		){
			argcheck_range("ignore",ignore,0.0,1.0);
			argcheck_gt("blocks",nblocks,(unsigned)1);

			double sumx=0.0;
			double sumx2=0.0;
			unsigned b;
			for(b=0;b<nblocks; b++){ // calculate average over each block
				double bav = getAverage(weight,
					Maths::min(1.0,ignore + double(b)*(1.0 - ignore)/double(nblocks) ) ,
					Maths::min(1.0,ignore + double(b+1)*(1.0 - ignore)/double(nblocks) ) );
				sumx += bav;
				sumx2 += Maths::sqr(bav);
			}
			return sqrt( (sumx2 - Maths::sqr(sumx)/double(nblocks))/(nblocks-1) );
	}


	Maths::Histogram1D_double MonitorBase::getHistogram(
		const MonitorBase &weight,
		double binwidth,
		unsigned start,
		unsigned end
		)const 
	{
		Maths::Histogram1D_double hist(binwidth);


		double sum = 0;
		double weightsum = 0;
		for(unsigned i = start;i<end;i++){
			sum += dataarray[i]*weight.getData(i);
			weightsum += weight.getData(i);
		}

		for(unsigned i = start;i<end;i++){
			hist.addPoint(dataarray[i],weight.getData(i));
		}
		return hist;
	}


	Maths::Histogram1D_double MonitorBase::getHistogram(
		const MonitorBase &weight,
		double binwidth, 
		unsigned start
		)const 
	{
		return getHistogram(weight, binwidth, start, dataarray.size());
	}

	void MonitorBase::setGranularity(unsigned newgran)
	{
		double order = ((double)newgran / (double)granularity);
		unsigned order_int = unsigned(order);
		if( (newgran <= 0) ||
			(!Maths::SigFigEquality<double>(order, order_int, 5 )) ){
				throw(ArgumentException("setGranularity(unsigned newgran): newgran must be a positive integer multiple of current granularity "));
		}

		if(nData() == 0 ){ // no data yet so just set granularity
			granularity *= order_int;
		}else{
			std::vector<double> newdataarray;
			gran_count = 0;
			gran_buffer = 0;
			for(size_t d_index = 0; d_index < dataarray.size(); d_index++ ){
				gran_buffer += dataarray[d_index]; // granuralise data by bundling it by the
				gran_count++; // ratio of new to old granularity (order_int)
				if(gran_count >= order_int){
					newdataarray.push_back(gran_buffer/(double)order_int);
					gran_buffer=0;
					gran_count=0;
				}
			}
			gran_buffer=0;
			gran_count=0;
			dataarray = newdataarray;
			granularity *= order_int;
		}
	}









	///////////////////////////////////////
	///
	/// MonitorContainer - Holds & Organises Monitors



	MonitorContainer::MonitorContainer()
	{
		name = "monc";
		reset();
	}

	int MonitorContainer::measure()
	{
		for(unsigned i=0;i<size();i++)
		{
			element(i).measure();
		}
		return 0;
	};

	int MonitorContainer::repeat()
	{
		for(unsigned i=0;i<size();i++)
		{
			element(i).repeat();
		}
		return 0;
	};

	void MonitorContainer::reset()
	{
		for(unsigned i=0;i<size();i++)
		{
			element(i).reset();
		}
	};

	void MonitorContainer::printHeader() const 
	{
		for(unsigned i=0;i<size();i++)
		{
			element(i).printHeader();
		}
	};

	void MonitorContainer::printCurData() const
	{
		for(unsigned i=0;i<size();i++)
		{
			element(i).printCurData();
		}
	}
	void MonitorContainer::printAllData() const
	{
		for(unsigned i=0;i<size();i++)
		{
			element(i).printAllData();
		}
	}
	void MonitorContainer::printAllDataParallel() const
	{
		if(size()==0){
			printf("No Monitors loaded\n");
			return;
		}
		printf("%s ",name.c_str());
		for(unsigned i=0;i<size();i++)
		{
			printf(" %10s",element(i).name.c_str());
		}
		printf("\n");
		for(unsigned idat=0;idat<element(0).nData();idat++)
		{
			printf("%s ",name.c_str());
			for(unsigned i=0;i<size();i++){
				element(i).printData(idat);
			}
			printf("\n");
		}
	}
}

