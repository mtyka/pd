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

#ifndef __HISTOGRAM_H
#define __HISTOGRAM_H

#include <deque>
#include <string>

namespace Maths
{
	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	///
	/// \author Mike Tyka  
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class Statistics{
	public:

		Statistics(){clear();}
		~Statistics(){};

		void clear(){
			m_n = 0;
			m_xlowest = 0; 
			m_xhighest = 0;
			m_xsum = 0;
			m_xsumsquared = 0;
		}

		void addPoint(double x){
			m_xsum += x;
			m_xsumsquared += sqr(x);

			if(m_n == 0) {
				m_xlowest = x;
				m_xhighest = x;
			} else {
				if(x < m_xlowest)
					m_xlowest = x;
				if(x > m_xhighest)
					m_xhighest = x;
			}
			m_n += 1;
		}
		double getAv() const{
			return (m_xsum / m_n);
		}
		double getVar() const{
			return fabs(((m_n * m_xsumsquared - sqr(m_xsum)) / sqr(m_n)));
		}
		double getStdDev()      const { return sqrt(getVar()); }
		double getHighest()     const { return m_xhighest;    }
		double getLowest()      const { return m_xlowest;     }
		double n()              const { return m_n;           }
		double getSum()         const { return m_xsum;        }
		double getSumSquared()	const { return m_xsumsquared; }

	protected:
		double m_n;           ///< Number of elements - can also be nonintegral
		double m_xlowest;     ///< Lowest Value encountered
		double m_xhighest;    ///< Largest Value encountered
		double m_xsum;        ///< Sum of m_data
		double m_xsumsquared; ///< Sum of squares
	};








	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	///
	/// \author Mike Tyka  
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	template <class T>
	class Histogram1D_Temp{
	public:

		Histogram1D_Temp();
		Histogram1D_Temp(double _xgridsize);
		Histogram1D_Temp(double X1, double X2, int Xsize);
		~Histogram1D_Temp();

		int init(double _xgridsize);
		int init(double X1, double X2, int Xsize);
		void clear();
		void print_row() const;
		void print_prob() const;
		void print_prob_cumul() const;
		void print_count() const;
		void print_count_double() const;
		void print_rdf() const;
		void print_gibbs(double Temp, double umk = 0.0) const;
		void print_gibbs_all(double Temp, double umk = 0.0) const;
		void add_point(double x){ add_point(x,T(1)); }
		void add_point(double x, const T &unitdata);
		void addPoint(double x){ add_point(x); }
		void addWidePoint(double x, double spread){ 
			double val;
			for( val = (x - spread); val <= (x + spread); val += m_xgridsize ) addPoint( val );
			add_point(x); 
		}
		void addPoint(double x, const T &unitdata){ add_point(x,unitdata); }

		double getAv() const{
			return (m_xsum / (double) (m_n + m_n_outside));
		}
		double getVar() const{
			return (((double) (m_n + m_n_outside) * m_xsumsquared - sqr(m_xsum)) / sqr((double) (m_n + m_n_outside)));
		}
		double getStdDev() const{
			return sqrt(getVar());
		}
		double getHighest() const{
			return m_xhighest;
		}
		double getLowest() const{
			return m_xlowest;
		}

		void getProb(std::deque<double> &prob) const{
			for(size_t i=0;i<m_data.size();i++){
				prob.push_back( m_data[i] / (double)m_n);
			}
		}

	public:  // public m_data
		std::string name;


		// inspectors exposing private member variables read-only
		T data(size_t index)  const { return m_data[index]; }
		size_t data_size()    const { return m_data.size(); }
		T n()                 const { return m_n;           }
		double x1()           const { return m_x1;          }
		double x2()					  const { return m_x2;          }
		double xgridsize()		const { return m_xgridsize;   }
		double xlowest()			const { return m_xlowest;     }
		double xhighest()		  const { return m_xhighest;    }
		double xsum()         const { return m_xsum;        }
		double xsumsquared()	const { return m_xsumsquared; }

	protected:
		std::deque<T> m_data; ///< Templatation allows non-integral counts
		T      m_n;           ///< Number of elements - can also be nonintegral
		T      m_n_outside;   ///< Number of elements outside of histogram range

		double m_x1;          ///< Boundaries of m_data field
		double m_x2;          ///< Boundaries of m_data field
		double m_xgridsize;   ///< Width of bins
		double m_xlowest;     ///< Lowest Value encountered
		double m_xhighest;    ///< Largest Value encountered
		double m_xsum;        ///< Sum of m_data
		double m_xsumsquared; ///< Sum of squares

		bool   m_staticsize;  ///< Does Histogram auto-grow to either side ?
	};



	// In 99% of case we just want an ordinary histogram - one that counts integers;
	// NOTE: shouldn't this be unsigned int ?

	typedef Histogram1D_Temp<int> Histogram1D;
	typedef Histogram1D_Temp<double> Histogram1D_double;

	// template Function Definitions

	template <class T>
	Histogram1D_Temp<T>::Histogram1D_Temp():name(""){
		init(1.0);
	}

	template <class T>
	Histogram1D_Temp<T>::Histogram1D_Temp(double _xgridsize):name(""){
		init(_xgridsize);
	}

	template <class T>
	Histogram1D_Temp<T>::Histogram1D_Temp(double X1, double X2, int Xsize):name(""){
		m_staticsize = true;
		init(X1, X2, Xsize);
	}



	template <class T>
	int Histogram1D_Temp<T>::init(double _xgridsize){
		m_x1=0;
		m_x2=0;
		m_xgridsize=_xgridsize;
		m_staticsize = false;
		clear();
		return 0;
	}

	template <class T>
	int Histogram1D_Temp<T>::init(double X1, double X2, int Xsize){
		if(Xsize > 0) {
			m_x1 = X1;
			m_x2 = X2;
			for(int i=0;i<Xsize;i++) m_data.push_back(0);
			m_xgridsize = (m_x2 - m_x1) / ((double) (m_data.size()));

		}else{
		}

		clear();
		return 0;
	}

	template <class T>
	Histogram1D_Temp<T>::~Histogram1D_Temp(){
	}

	template <class T>
	void Histogram1D_Temp<T>::clear(){
		unsigned i;
		for(i=0;i<m_data.size();i++) m_data[i]=0;
		m_n = 0;
		m_n_outside = 0;
		m_xsum = 0;
		m_xsumsquared = 0;
	}

	template <class T>
	void Histogram1D_Temp<T>::print_row() const{
		int x;

		if(m_data.size() > 0) {

			for(x = 0; x < m_data.size(); x++) {
				printf("%6d\t", int(m_data[x]));
			}
			printf("\n\n");
			for(x = 0; x < m_data.size(); x++) {
				printf("%s %8.3lf\t", 
					name.c_str(),
					(double) m_data[x] / (double) m_n);
			}
			printf("\nTotal Points: %d (%d)\n", int(m_n), int(m_n_outside));
		}
	}

	//template <class T>
	//void Histogram1D_Temp<T>::print_vertical(FILE * file, const std::string &prefix, int greaterthen) const{
	//	fprintf(file, "%sav %e \n", prefix.c_str(), getAverageSubdata());
	//	//fprintf(file,"%svar %e \n",prefix.c_str(),variance());
	//	fprintf(file, "%sstd %e \n", prefix.c_str(), standarddev());

	//	if(m_data.size() <= 0)
	//		return;
	//	int x, start = 0, end = (int)m_data.size() - 1;
	//	for(start = 0; start < (int)m_data.size(); start++)
	//		if(m_data[start] >= greaterthen)
	//			break;
	//	for(end = (int)m_data.size() - 1; end > 0; end--)
	//		if(m_data[end] >= greaterthen)
	//			break;
	//	for(x = start; x <= end; x++) {

	//		fprintf(file, "%shst%10.5lf\t%10.8lf\t%8d\n",
	//			prefix.c_str(),double( m_x1 + (double) x * (double) m_xgridsize), double( int(m_data[x]) / (double)int(m_n)), int(m_data[x]));
	//	}
	//	fprintf(file, "%shsn %d %d)\n",prefix.c_str(), int(m_n), int(m_n_outside));
	//}

	template <class T>
	void Histogram1D_Temp<T>::print_prob() const{

		if(m_data.size() <= 0) return;
		int x, start = 0, end = (int)m_data.size() - 1;
		for(x = start; x <= end; x++) {
			printf("%s %8.7f\t%8.7lf\n",
				name.c_str(),
				m_x1 + (double) x * (double) m_xgridsize,
				(((double) m_data[x]) / (double) m_n));
		}
	}

	template <class T>
	void Histogram1D_Temp<T>::print_prob_cumul() const{
		double sum=0.0;
		if(m_data.size() <= 0) return;
		int x, start = 0, end = (int)m_data.size() - 1;
		for(x = start; x <= end; x++) {
			sum+=(double) m_data[x];
			printf("%s %8.7f\t%8.7lf\n",
				name.c_str(),
				m_x1 + (double) x * (double) m_xgridsize,
				((sum) / (double) m_n));
		}
	}


	template <class T>
	void Histogram1D_Temp<T>::print_count() const{

		if(m_data.size() <= 0) return;
		int x, start = 0, end = (int)m_data.size() - 1;
		for(x = start; x <= end; x++) {
			printf("%s %8.4f\t%10d\n",
				name.c_str(),
				m_x1 + (double) x * (double) m_xgridsize, 
				int(m_data[x]));
		}
	}

	template <class T>
	void Histogram1D_Temp<T>::print_count_double() const{

		if(m_data.size() <= 0) return;
		int x, start = 0, end = (int)m_data.size() - 1;
		for(x = start; x <= end; x++) {
			printf("%s %8.4f\t%10.7e\n",
				name.c_str(),
				m_x1 + (double) x * (double) m_xgridsize, 
				double(m_data[x]));
		}
	}

	template <class T>
	void Histogram1D_Temp<T>::print_rdf() const{

		if(m_data.size() <= 0) return;
		int x, start = 0, end = (int)m_data.size() - 1;
		for(x = start; x <= end; x++) {
			printf("%s %8.7f\t%8.7lf\n",
				name.c_str(),
				m_x1 + (double) x * (double) m_xgridsize,
				(((double) m_data[x]) / (double) m_n / (4*MathConst::PI*sqr(m_x1 + (double) x * (double) m_xgridsize)) ));
		}
	}

	template <class T>
	void Histogram1D_Temp<T>::print_gibbs(double Temp, double umk) const{

		if(m_data.size() <= 0) return;
		int x, start = 0, end = (int)m_data.size() - 1;
		for(x = start; x <= end; x++) {
			printf("%s %8.4f\t%8.3lf\t%8.3lf\n",
				name.c_str(),
				m_x1 + (double) x * (double) m_xgridsize,
				-0.0019792493288573747*Temp*log(((double) m_data[x]) / (double) m_n),
				-0.5*umk*(m_x1 + ((double) x + 0.5)* (double) m_xgridsize)
				- 0.0019792493288573747*Temp*log(((double) m_data[x]) / (double) m_n)
				);
		}
	}

	template <class T>
	void Histogram1D_Temp<T>::print_gibbs_all(double Temp, double umk) const{

		if(m_data.size() <= 0) return;
		int x, start = 0, end = (int)m_data.size() - 1;
		for(x = start; x <= end; x++) {
			printf("%s %8.4f\t%10d\t%8.3lf\t%8.3lf\t%8.3lf\n",
				name.c_str(),
				m_x1 + (double) x * (double) m_xgridsize,
				int(m_data[x]),
				(((double) m_data[x] + 0.5) / (double) m_n),
				-0.0019792493288573747*Temp*log(((double) m_data[x]) / (double) m_n),
				-0.5*umk*(m_x1 + ((double) x + 0.5)* (double) m_xgridsize)
				- 0.0019792493288573747*Temp*log(((double) m_data[x]) / (double) m_n)
				);
		}
	}

	template <class T>
	void Histogram1D_Temp<T>::add_point(double x, const T &unitdata){
		int xc;
		m_xsum += x;
		m_xsumsquared += sqr(x);

		if(m_n == 0) {
			m_xlowest = x;
			m_xhighest = x;
		} else {
			if(x < m_xlowest)
				m_xlowest = x;
			if(x > m_xhighest)
				m_xhighest = x;
		}

		if(m_data.size()==0){// are there any bins yet ?
			// create the first bin:

			double sigfig=1;
			//do{
			m_x1=floor(sigfig*x)/sigfig; // set beginning of first bin to the first value
			sigfig*=10;
			//}while((m_x1+m_xgridsize)<x);
			m_x2=m_x1+m_xgridsize;
			m_data.push_back(0);
			add_point(x,unitdata);
			return;
		}else{

			if((x < m_x1) || (x >= m_x2)) {

				if(m_staticsize){
					m_n_outside++;
				}else{
					//printf("ext size \n");
					int newbins;
					if(x<m_x1){
						newbins = int(ceil((m_x1-x)/m_xgridsize) + 1);
						// runtime assertion to prevent excessive memory usage explosion
						ASSERT(newbins * sizeof(T) < 1000000000,ArgumentException,"Runtime ASSERT failed - the histogram you are trying create is too fine grained\nfor the scale of m_data entered - check the m_data is ok or increase the binwidth");
						for(int i=0;i<newbins;i++) m_data.push_front(0);
						m_x1 -= newbins*m_xgridsize;
					}else{
						newbins = int(ceil((x-m_x2)/m_xgridsize) + 1);
						ASSERT(newbins * sizeof(T) < 1000000000,ArgumentException,"Runtime ASSERT failed - the histogram you are trying create is too fine grained\nfor the scale of m_data entered - check the m_data is ok or increase the binwidth");
						for(int i=0;i<newbins;i++) m_data.push_back(0);
						m_x2 += newbins*m_xgridsize;
					}
					add_point(x,unitdata); // try again to add point
				}
				return;

			}
			x -= m_x1;
			x /= (m_xgridsize);
			xc = (int) floor(x);
			if(xc < 0) return;
			m_data[xc] += unitdata;
		}

		m_n += unitdata;
	}





















	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	///
	/// \author Mike Tyka  
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	template <class T>
	class Histogram2D_Temp{
	public:
		std::deque< std::deque< T > > m_data;
		int m_n, m_n_outside;

		double m_x1, m_x2;
		double m_xgridsize;
		double m_xlowest;
		double m_xhighest;
		double m_xsum;
		double m_xsumsquared;

		double y1, y2;
		double ygridsize;
		double ylowest;
		double yhighest;
		double ysum;
		double ysumsquared;

		bool m_staticsize;

	public:

		Histogram2D_Temp();
		Histogram2D_Temp(double _xgridsize,double _ygridsize);
		Histogram2D_Temp(const Histogram2D_Temp &copyhist);
		Histogram2D_Temp(double X1, double X2, int Xsize,
			double Y1, double Y2, int Ysize);
		~Histogram2D_Temp();

#ifndef SWIG
		Histogram2D_Temp& operator = (const Histogram2D_Temp &copyhist);
#endif

		int init(double _xgridsize,double _ygridsize);
		int init(double X1, double X2, int Xsize, double Y1, double Y2, int Ysize);
		void clear();
		void print_map_prob();
		void print_map_gibbs(double Temp, bool raw = false);
		void print_vertical(FILE * file, const std::string &prefix);
		void add_point(double x, double y){ add_point(x,y,T(1)); }
		void add_point(double x, double y, const T &unitdata);

		double getAverageSubdata(){
			return (m_xsum / (double) (m_n + m_n_outside));
		}
		double variance(){
			return (((double) (m_n + m_n_outside) * m_xsumsquared - sqr(m_xsum)) / sqr((double) (m_n + m_n_outside)));
		}
		double standarddev(){
			return sqrt(variance());
		}
		double highest(){
			return m_xhighest;
		}
		double lowest(){
			return m_xlowest;
		}
	};



	/// In 99% of case we just want an ordinary histogram - one that counts integers;
	/// NOTE: shouldn't this be unsigned int ?

	typedef Histogram2D_Temp<int> Histogram2D;

	// template Function Definitions

	template <class T>
	Histogram2D_Temp<T>::Histogram2D_Temp(){
		init(1.0,1.0);
	}

	template <class T>
	Histogram2D_Temp<T>::Histogram2D_Temp(double _xgridsize, double _ygridsize){
		init(_xgridsize, _ygridsize);
	}

	template <class T>
	Histogram2D_Temp<T>::Histogram2D_Temp(double X1, double X2, int Xsize,
		double Y1, double Y2, int Ysize){
			m_staticsize = true;
			init(X1, X2, Xsize, Y1, Y2, Ysize);
	}

	template <class T>
	Histogram2D_Temp<T>::Histogram2D_Temp(const Histogram2D_Temp &copyhist){
		m_staticsize = true;
		(*this) = copyhist;
	}

	template <class T>
	Histogram2D_Temp<T>& Histogram2D_Temp<T>::operator = (const Histogram2D_Temp &copyhist){
		m_n = copyhist.m_n;
		m_n_outside = copyhist.m_n_outside;
		m_x1 = copyhist.m_x1;
		m_x2 = copyhist.m_x2;
		m_xgridsize = copyhist.m_xgridsize;
		m_xlowest = copyhist.m_xlowest;
		m_xhighest = copyhist.m_xhighest;
		m_xsum = copyhist.m_xsum;
		m_xsumsquared = copyhist.m_xsumsquared;
		y1 = copyhist.y1;
		y2 = copyhist.y2;
		ygridsize = copyhist.ygridsize;
		ylowest = copyhist.ylowest;
		yhighest = copyhist.yhighest;
		ysum = copyhist.ysum;
		ysumsquared = copyhist.ysumsquared;
		m_data = copyhist.m_data;
		m_staticsize = copyhist.m_staticsize;
		return (*this);
	}

	template <class T>
	int Histogram2D_Temp<T>::init(double _xgridsize, double _ygridsize){
		m_x1=0;
		m_x2=0;
		y1=0;
		y2=0;
		m_xgridsize=_xgridsize;
		ygridsize=_ygridsize;
		m_staticsize = false;
		clear();
		return 0;
	}

	template <class T>
	int Histogram2D_Temp<T>::init(double X1, double X2, int Xsize,
		double Y1, double Y2, int Ysize){
			if((Xsize > 0)&&(Ysize > 0)) {
				m_x1 = X1;
				m_x2 = X2;
				y1 = Y1;
				y2 = Y2;
				for(int i=0;i<Xsize;i++){
					m_data.push_back( std::deque < T > () );
					for(int j=0;j<Ysize;j++){
						m_data[i].push_back(0);
					}
				}
				m_xgridsize = (m_x2 - m_x1) / ((double) (m_data.size()));
				ygridsize = (y2 - y1) / ((double) (m_data[0].size()));
			}else{
			}

			clear();
			return 0;
	}

	template <class T>
	Histogram2D_Temp<T>::~Histogram2D_Temp(){
	}

	template <class T>
	void Histogram2D_Temp<T>::clear(){
		unsigned i,j;
		for(i=0;i<m_data.size();i++){
			for(j=0;j<m_data[i].size();j++){
				m_data[i][j]=0;
			}
		}
		m_n = 0;
		m_n_outside = 0;
		m_xsum = 0;
		m_xsumsquared = 0;
	}


	/// prints a map of probabilities i.e. dividing the histogram count by the total number of points
	/// inside the histogram
	template <class T>
	void Histogram2D_Temp<T>::print_map_prob(){
		int x,y;

		if(m_data.size() > 0) {

			printf(" \t");
			for(y = 0; y < m_data[0].size(); y++) {
				printf("%8.3lf\t", y1 + y*ygridsize );
			}
			printf("\n \t");
			for(y = 0; y < m_data[0].size(); y++) {
				printf("--------\t", y1 + y*ygridsize );
			}
			printf("\n");
			for(x = 0; x < m_data.size(); x++) {
				printf("%8.3lf |\t", m_x1 + x*m_xgridsize );
				for(y = 0; y < m_data[x].size(); y++) {
					printf("%8.3lf\t", (double) m_data[x][y] / (double) m_n);
				}
				printf("\n");
			}
			printf("\nTotal Points: %d (%d)\n", int(m_n), int(m_n_outside));
		}else{
			printf("--empty--");
		}
	}

	/// prints map converting probabilities into free energies using A = -RT ln (P)
	/// where T is set by the parameter 'Temp' (in kcal)
	template <class T>
	void Histogram2D_Temp<T>::print_map_gibbs(double Temp, bool raw){
		int x,y;

		if(m_data.size() > 0) {
			if(!raw){
				printf(" \t");
				for(y = 0; y < m_data[0].size(); y++) {
					printf("%8.3lf\t", y1 + y*ygridsize );
				}
				printf("\n \t");
				for(y = 0; y < m_data[0].size(); y++) {
					printf("--------\t", y1 + y*ygridsize );
				}
				printf("\n");
			}
			for(x = 0; x < m_data.size(); x++) {
				if(!raw) printf("%8.3lf |\t", m_x1 + x*m_xgridsize );
				for(y = 0; y < m_data[x].size(); y++) {
					printf("%8.3lf\t", -0.0019792493288573747*Temp*log(((double) m_data[x][y] + 0.5) / (double) m_n));
				}
				printf("\n");
			}
			//printf("\nTotal Points: %d (%d)\n", m_n, m_n_outside);
		}else{
			if(!raw){
				printf("--empty--");
			}
		}
	}

	template <class T>
	void Histogram2D_Temp<T>::print_vertical(FILE * file, const std::string &prefix){
		fprintf(file, "%sav %e \n", prefix.c_str(), getAverageSubdata());
		//fprintf(file,"%svar %e \n",prefix.c_str(),variance());
		fprintf(file, "%sstd %e \n", prefix.c_str(), standarddev());

		if(m_data.size() <= 0)
			return;
		size_t x;
		for(x = 0; x < m_data.size(); x++) {
			for(size_t y = 0; y < m_data[x].size(); y++) {

				fprintf(file, "%shst%10.5lf\t%10.5lf\t%10.8lf\t%8d\n",
					prefix.c_str(),
					m_x1 + (double) x * (double) m_xgridsize,
					y1 + (double) y * (double) ygridsize,
					(double) int(m_data[x][y]) / (double) m_n, int(m_data[x][y]));
			}
		}
		fprintf(file, "%shsn %d %d)\n",prefix.c_str(), int(m_n), int(m_n_outside));
	}

	template <class T>
	void Histogram2D_Temp<T>::add_point(double x, double y, const T &unitdata){
		int xc;
		int yc;
		m_xsum += x;
		m_xsumsquared += sqr(x);

		if(m_n == 0) {
			m_xlowest = x;
			m_xhighest = x;
			ylowest = y;
			yhighest = y;
		} else {
			if(x < m_xlowest) m_xlowest = x;
			if(x > m_xhighest) m_xhighest = x;
			if(y < ylowest) ylowest = y;
			if(y > yhighest) yhighest = y;
		}

		if(m_data.size()==0){// are there any bins yet ?
			// create the first bin:

			double sigfig=1;
			m_x1=floor(sigfig*x)/sigfig; // set beginning of first bin to the first value
			sigfig*=10;
			m_x2=m_x1+m_xgridsize;
			y1=floor(sigfig*y)/sigfig; // set beginning of first bin to the first value
			sigfig*=10;
			y2=y1+ygridsize;
			m_data.push_back( std::deque < T > () );
			m_data[0].push_back( 0 );
			add_point(x,y,unitdata);
			return;
		}else{

			if((x < m_x1) || (x >= m_x2) ||
				(y < y1) || (y >= y2)) {

					if(m_staticsize){
						m_n_outside++;
					}else{
						//printf("ext size \n");
						int newbins;
						if(x<m_x1){
							newbins = int(ceil((m_x1-x)/m_xgridsize) + 1);
							for(int i=0;i<newbins;i++){
								std::deque < T > newcol;
								for(int j=0;j<m_data[0].size();j++){
									newcol.push_back( 0 );
								}
								m_data.push_front(newcol);
							}
							m_x1 -= newbins*m_xgridsize;
						}else
							if(x>m_x2){
								newbins = int(ceil((x-m_x2)/m_xgridsize) + 1);
								for(int i=0;i<newbins;i++){
									std::deque < T > newcol;
									for(int j=0;j<m_data[0].size();j++){
										newcol.push_back( 0 );
									}
									m_data.push_back(newcol);
								}
								m_x2 += newbins*m_xgridsize;
							}else
								if(y<y1){
									newbins = int(ceil((y1-y)/ygridsize) + 1);
									for(int j=0;j<m_data.size();j++){
										for(int i=0;i<newbins;i++) m_data[j].push_front(0);
									}
									y1 -= newbins*ygridsize;
								}else
									if(y>y2){
										newbins = int(ceil((y-y2)/ygridsize) + 1);
										for(int j=0;j<m_data.size();j++){
											for(int i=0;i<newbins;i++) m_data[j].push_back(0);
										}
										y2 += newbins*ygridsize;
									}


									add_point(x,y,unitdata); // try again to add point
					}
					return;

			}
			x -= m_x1;
			x /= (m_xgridsize);
			xc = (int) floor(x);
			if(xc < 0) return;
			y -= y1;
			y /= (ygridsize);
			yc = (int) floor(y);
			if(yc < 0) return;
			m_data[xc][yc] += unitdata;
		}

		m_n++;
	}
}

#ifdef SWIG
%template(Histogram1D) Maths::Histogram1D_Temp<int>;
%template(Histogram1D_double) Maths::Histogram1D_Temp<double>;
%template(Histogram2D) Maths::Histogram2D_Temp<int>;
#endif

#endif

