#include <iostream>
#include <process.h>

#include <process.h>
#include <boost/random.hpp> 

#include "hmm.h"
#include "functions.h"

int main() {
	    using namespace boost::numeric::ublas;
    matrix<double> m (3, 4);
    for (unsigned i = 0; i < m.size1 (); ++ i)
        for (unsigned j = 0; j < m.size2 (); ++ j)
            m (i, j) = 4 * i + j;
    std::cout << m << std::endl;

		getchar();
		//return 0;


  std::string path = "C:\\Users\\Estudiante\\Documents\\GitHub\\msc-thesis-code\\";
  // std::string path = "E:\\ESCUELA\\CIMAT\\4 Semestre\\ST2\\prog\\"
	
  auto data = readCSV<uint>(path + "pruebas\\cuervo1f_120.csv");

	//data.resize(4, true);

  uint nn = 6;
  uint kk = 120;    
  uint tt = data.size();

  HMM::params p(nn, kk, true);

  HMM h(data);
  h.EM(p);

  auto res = p.hidden;
  writeCSV<uint>(res, path + "pruebas\\cuervo1f_mground.csv");

  getchar();
  return 0;
}