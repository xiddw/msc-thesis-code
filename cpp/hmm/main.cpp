#include <iostream>
#include <process.h>

#include <process.h>
#include <boost/random.hpp> 

#include "hmm.h"
#include "functions.h"

int main() {

  HMM hmm(20, 15);
  hmm.genKeys(6);

  uint nn = 3;
  uint kk = 5;  
  
  HMM::params p = HMM::generateParams(nn, kk);
  std::cout << std::endl << "P: " << p << std::endl;

  HMM::params q(nn, kk);
  std::cout << std::endl << "Q: " << q << std::endl;

  HMM::params r(nn, kk, true);
  std::cout << std::endl << "R: " << r << std::endl;

  std::cout << "Hola" << std::endl;
  auto d = hmm.sample(r, 100);
  std::cout << d << std::endl;
  
  getchar();
  return 0;
}