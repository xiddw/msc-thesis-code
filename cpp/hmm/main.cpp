#include <iostream>
#include <process.h>

#include <process.h>
#include <boost/random.hpp> 

#include <ctime>
#include <omp.h>

#include "hmm.h"
#include "functions.h"



int main() {
  // std::string path = "C:\\Users\\Estudiante\\Documents\\GitHub\\msc-thesis-code\\";
  std::string path = "E:\\ESCUELA\\CIMAT\\4 Semestre\\ST2\\prog\\";
	
  auto data = readCSV<uint>(path + "pruebas\\cuervo1f_120.csv");

  uint nn = 6;
  uint kk = 120;    
  uint tt = data.size();

  HMM::params m(nn, kk);

  HMM h(data);
  
	int max = 2*h.MaxIterHMM();

	int i, tid, nthreads;

	// Auxiliar variables for omp process :)
	double *ll = new double[max], 
					mm1, mm2;

	HMM::params *pp = new HMM::params[max], 
							 oo1, oo2;

	clock_t ts = clock();

	#pragma omp parallel for private(tid, i) shared(ll, pp) 
	for(i=0; i<max; ++i) {
		bool s = (i >= max/2);
		HMM::params p(nn+s, kk, true);
		ll[i] = h.EM(p);
		pp[i] = p;

		std::cout << "."; 
	}

	mm1 = -DBL_MAX;
	mm2 = -DBL_MAX;

	int ii1 = 0, ii2;
	for(i=0; i<max/2; ++i) {
		if(ll[i] > mm1) {
			mm1 = ll[i];
			ii1 = i;
		}
	}
	for(; i<max; ++i) {
		if(ll[i] > mm2) {
			mm2 = ll[i];
			ii2 = i;
		}
	}

	oo1 = pp[ii1];
	oo2 = pp[ii2];

	clock_t te = clock();

	printf("\n Time taken: %.2fs\n", (double)(te -ts)/CLOCKS_PER_SEC);

  std::cout << "best overal [model 1] (iter=>" << ii1 << "): " << mm1 << std::endl;
	std::cout << "best overal [model 2] (iter=>" << ii2 << "): " << mm2 << std::endl;

  auto res1 = oo1.hidden;
	auto res2 = oo2.hidden;
  writeCSV<uint>(res1, path + "pruebas\\cuervo1f_" + std::to_string(nn) + "_mground.csv");
	writeCSV<uint>(res2, path + "pruebas\\cuervo1f_" + std::to_string(nn+1) + "_mground.csv");

  getchar();
  return 0;
}